#!/usr/bin/env python2
"""
Retrieves X, Y, Z components from SEIS files.
Saves to HDF5 file for all grid points.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 2 June 2016

ISSUES: All options of EMOD3D are not considered. eg: scale.
"""

from glob import glob
from math import sin, cos, radians
from struct import unpack

from mpi4py import MPI
import numpy as np

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank() # process number
size = comm.Get_size() # total process number

# H5 parallel IO requires every process to create every data structure
# it makes no sense and too much additional work
if rank == 0:
    from os import stat
    from sys import byteorder
    import h5py as h5
    h5_file = 'virtual.hdf5'
    h5p = h5.File(h5_file, 'w')
    # check for messages in queue
    peek = comm.Iprobe
    ANY_SOURCE = MPI.ANY_SOURCE

# 8 byte c-string has up to 7 characters
STAT_CHAR = 8
# a station with a name of this length is virtual
VSTAT_LEN = 4
# number of bytes per value, int/float = 4 Byte on 64 bit machines
SIZE_INT = 4
SIZE_FLT = 4
SIZE_SEISHEAD = SIZE_INT * 5 + SIZE_FLT * 5 + STAT_CHAR
# binary format properties
N_COMPS = 9
# values of interest in components, index in N_COMPS, description
# changing these requires changing logic
MY_COMPS = {0:'090', 1:'000', 2:'VUP'}
N_MY_COMPS = len(MY_COMPS)
# MPI TAGS
TAG_NSTAT = 73
TAG_STATI = 85
TAG_STATC = 182

# process 0 will stop and check for/process messages
def recv_interrupt():
    while peek(ANY_SOURCE, TAG_NSTAT):
        num_stat = comm.recv(source = ANY_SOURCE, tag = TAG_NSTAT)
        try:
            h5p.attrs['NSTAT'] += num_stat
        except KeyError:
            h5p.attrs['NSTAT'] = num_stat
        print 'running total nstat = %d' % h5p.attrs['NSTAT']

    while peek(ANY_SOURCE, TAG_STATI):
        name, x, y, lat, lon = comm.recv(source = ANY_SOURCE, tag = TAG_STATI)
        try:
            dset = h5p.create_dataset(name)
            dset.attrs['X'] = x
            dset.attrs['Y'] = y
            dset.attrs['LAT'] = lat
            dset.attrs['LON'] = lon
        except:
            pass

    while peek(ANY_SOURCE, TAG_STATC):
        data = np.zeros(20000 * 3, shape = (20000, 3))
        comm.Recv(data, source = ANY_SOURCE, tag = TAG_STATC)

# automatic endianness detection
# first assume that endianness is native, check if filesize matches
def get_swap(fp):
    # number of stations in seis file
    fp.seek(0)
    ns = unpack('i', fp.read(SIZE_INT))[0]
    # first station NT value (should all be same)
    fp.seek(SIZE_INT + 4 * SIZE_INT)
    nt = unpack('i', fp.read(SIZE_INT))[0]

    # assuming correct values read, this is the filesize
    fs = SIZE_INT + ns * (SIZE_SEISHEAD + SIZE_FLT * N_COMPS * nt)
    
    if fs == stat(fp.name).st_size:
        return False
    return True

seis_file_list = glob('OutBin/*_seis-?????.e3d')
rank_file_list = seis_file_list[rank::size]

# byte swap determination
if rank == 0:
    fp = open(seis_file_list[0], 'r')
    if get_swap(fp):
        if byteorder == 'little':
            # big endian src, little endian local
            int_s = '>i'
            flt_s = '>f'
        else:
            # little endian src, big endian local
            int_s = '<i'
            flt_s = '<f'
    else:
        # no byte swapping
        int_s = 'i'
        flt_s = 'f'
    fp.close()
    isfs = [int_s, flt_s]
else:
    isfs = None
INT_S, FLT_S = comm.bcast(isfs, root = 0)

# common parameters
if rank == 0:
    fp = open(seis_file_list[0], 'rb')
    # read what is assumed to be common amongst stations from first station
    fp.seek(SIZE_INT * 5)
    nt = unpack(INT_S, fp.read(SIZE_INT))[0]
    dt = unpack(FLT_S, fp.read(SIZE_FLT))[0]
    hh = unpack(FLT_S, fp.read(SIZE_FLT))[0]
    rot = unpack(FLT_S, fp.read(SIZE_FLT))[0]
    fp.close()
    # also save to HDF file
    h5p.attrs['NSTAT'] = 0
    h5p.attrs['NT'] = nt
    h5p.attrs['DT'] = dt
    h5p.attrs['HH'] = hh
    h5p.attrs['ROT'] = rot
    # pass to all other processes
    specs = (nt, dt, hh, rot)
else:
    specs = None
nt, dt, hh, rot = comm.bcast(specs, root = 0)

# rotation matrix for converting to 090, 000
# ver is inverted (* -1)
theta = radians(rot)
rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
        [-sin(theta), -cos(theta), 0], \
        [0, 0, -1]])

for si, seis_file in enumerate(rank_file_list):
    print('[R%d] processing %d of %d...' \
            % (rank, rank + si * size + 1, len(seis_file_list)))
    fp = open(seis_file, 'rb')
    # speed optimise
    seek = fp.seek
    read = fp.read
    # shortcuts for processing input
    read_int = lambda : unpack(INT_S, read(SIZE_INT))[0]
    read_flt = lambda : unpack(FLT_S, read(SIZE_FLT))[0]

    seek(0)
    num_stat = read_int()
    comm.send(num_stat, dest = 0, tag = TAG_NSTAT)
    if rank == 0:
        recv_interrupt()

    # read at once, all component data below header as float array
    # major speedup compared to manual processing
    seek(SIZE_INT + num_stat * SIZE_SEISHEAD)
    comp_data = np.fromfile(fp, dtype = FLT_S)
    if rank == 0:
        recv_interrupt()

    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        if len(str(read(STAT_CHAR)).rstrip('\0')) == VSTAT_LEN:
            # station is a virtual entry

            seek( - STAT_CHAR - 5 * SIZE_FLT - 4 * SIZE_INT, 1)
            x = read_int()
            y = read_int()

            seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            stat_params = ('%s%s' % (str(x).zfill(4), str(y).zfill(4)), \
                    x, y, read_flt(), read_flt())
            comm.send(stat_params, dest = 0, tag = TAG_STATI)

            comps = np.dot(np.dstack(( \
                    comp_data[0 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[1 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[2 + stat_i * N_COMPS : : num_stat * N_COMPS])), rot_matrix)
            comm.Send([comps, MPI.FLOAT], dest = 0, tag = TAG_STATC)
            if rank == 0:
                recv_interrupt()

    fp.close()


print('[R%s] processes complete.' % rank)
comm.barrier()
if rank == 0:
    recv_interrupt()
    h5p.close()
print('[R%s] exit.' % rank)

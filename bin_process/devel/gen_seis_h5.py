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
from os import stat, remove, path
from struct import unpack
from sys import byteorder

from mpi4py import MPI
import numpy as np
import h5py as h5

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

h5_file = 'virtual.hdf5'
if path.isfile(h5_file):
    remove(h5_file)
# HDF5 throws errors when opened in 'w' mode
# happens when creating some groups (not duplicates)
h5p = h5.File(h5_file, 'a')
seis_file_list = glob('OutBin/*_seis-?????.e3d')

fp = open(seis_file_list[0], 'rb')
if get_swap(fp):
    if byteorder == 'little':
        # big endian src, little endian local
        INT_S = '>i'
        FLT_S = '>f'
    else:
        # little endian src, big endian local
        INT_S = '<i'
        FLT_S = '<f'
else:
    # no byte swapping
    INT_S = 'i'
    FLT_S = 'f'
fp.close()

for si, seis_file in enumerate(seis_file_list):
    print('processing %d of %d...' \
            % (si + 1, len(seis_file_list)))
    fp = open(seis_file, 'rb')
    # speed optimise
    seek = fp.seek
    read = fp.read

    # shortcuts for processing input
    read_int = lambda : unpack(INT_S, read(SIZE_INT))[0]
    read_flt = lambda : unpack(FLT_S, read(SIZE_FLT))[0]

    num_stat = read_int()

    # write common data to HDF5
    try:
        h5p.attrs['NSTAT'] += num_stat
    except KeyError:
        h5p.attrs['NSTAT'] = num_stat
        # read what is assumed to be common amongst stations from first station
        seek(SIZE_INT * 5)
        nt = read_int()
        dt = read_flt()
        hh = read_flt()
        rot = read_flt()
        # rotation matrix for converting to 090, 000
        # ver is inverted (* -1)
        theta = radians(rot)
        rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
                [-sin(theta), -cos(theta), 0], \
                [0, 0, -1]])
        h5p.attrs['NT'] = nt
        h5p.attrs['DT'] = dt
        h5p.attrs['HH'] = hh
        h5p.attrs['ROT'] = rot
        h5p.attrs['NCOMPS'] = N_MY_COMPS
        for n, key in enumerate(MY_COMPS):
            h5p.attrs['COMP_%d' % (n)] = MY_COMPS[key]


    # read at once, all component data below header as float array
    # major speedup compared to manual processing
    seek(SIZE_INT + num_stat * SIZE_SEISHEAD)
    comp_data = np.fromfile(fp, dtype = FLT_S)

    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        stat = str(read(STAT_CHAR)).rstrip('\0')
        if len(stat) == VSTAT_LEN:
            # station is a virtual entry

            seek( - STAT_CHAR - 5 * SIZE_FLT - 4 * SIZE_INT, 1)
            x = read_int()
            y = read_int()

            # station group in HDF5 (XXXXYYYY)
            g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))
            try:
                s_group = h5p.create_group(g_name)
            except ValueError, e:
                # if group created previously
                print 'duplicate station:', g_name
                #continue
                s_group = h5p[g_name]

            s_group.attrs['NAME'] = stat
            s_group.attrs['X'] = x
            s_group.attrs['Y'] = y
            seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            s_group.attrs['LAT'] = read_flt()
            s_group.attrs['LON'] = read_flt()

            comps = np.dot(np.dstack(( \
                    comp_data[0 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[1 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[2 + stat_i * N_COMPS : : num_stat * N_COMPS])), rot_matrix)

            try:
                h5p.create_dataset(g_name + '/VEL', (nt, N_MY_COMPS), \
                        dtype='f', data = comps)
            except RuntimeError:
                # this station had beed added before
                h5p[g_name + '/VEL'][...] = comps

    fp.close()

h5p.close()

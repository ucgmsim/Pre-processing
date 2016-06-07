#!/usr/bin/env python2
"""
Retrieves X, Y, Z components from SEIS files.
Saves to HDF5 file for all grid points.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 2 June 2016

ISSUES: All options of EMOD3D are not considered. eg: scale.
"""

from array import array
from itertools import izip
from math import sin, cos, radians
from os import stat
from struct import unpack
from sys import byteorder

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

seis_file = 'LPSim-2011Feb22b560_v1_Cantv1_64-h0.100_v3.04_seis-01669.e3d'
h5_file = 'virtual.hdf5'
h5p = h5.File(h5_file, 'a')

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
    
    if fs == stat(seis_file).st_size:
        return False
    return True

fp = open(seis_file, 'rb')
# whether in non-native endian
swap = get_swap(fp)
# speed optimise
seek = fp.seek
read = fp.read

# python format string, can handle byte swapping
if swap:
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

# shortcuts for processing input
read_int = lambda : unpack(INT_S, read(SIZE_INT))[0]
read_flt = lambda : unpack(FLT_S, read(SIZE_FLT))[0]

seek(0)
num_stat = read_int()
# read what is assumed to be common amongst stations from first station
seek(SIZE_INT * 5)
nt = read_int()
dt = read_flt()
hh = read_flt()
rot = read_flt()

# write common data to HDF5
try:
    h5p.attrs['NSTAT'] += num_stat
except KeyError:
    h5p.attrs['NSTAT'] = num_stat
    h5p.attrs['NT'] = nt
    h5p.attrs['DT'] = dt
    h5p.attrs['HH'] = hh
    h5p.attrs['ROT'] = rot

# read at once, all component data below header as float array
# major speedup compared to manual processing
seek(SIZE_INT + num_stat * SIZE_SEISHEAD)
# below method would be faster if could disregard the rest
#dt = np.dtype([ 
#        ('c1', FLT_S),
#        ('c2', FLT_S),
#        ('c3', FLT_S),
#        ('rest', np.void(SIZE_FLT * 6))])
comp_data = np.fromfile(fp, dtype = FLT_S)
# rotation matrix for converting to 000, 090
# ver is inverted (* -1)
theta = radians(rot)
rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
        [-sin(theta), -cos(theta), 0], \
        [0, 0, -1]])

for stat_i in xrange(num_stat):
    seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
    if len(str(read(STAT_CHAR)).rstrip('\0')) == VSTAT_LEN:
        # station is a virtual entry

        seek( - STAT_CHAR - 5 * SIZE_FLT - 4 * SIZE_INT, 1)
        x = read_int()
        y = read_int()
        seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)

        # station group in HDF5 (XXXXYYYY)
        g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))
        try:
            s_group = h5p.create_group(g_name)
        except ValueError:
            # most likely group created previously
            s_group = h5p[g_name]

        s_group.attrs['X'] = x
        s_group.attrs['Y'] = y
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

#!/usr/bin/env python2
"""
Retrieves X, Y, Z components from SEIS files.
Processes ready for BB.
Saves to Binary, ready for BB.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 13 September 2016

ISSUES: All options of EMOD3D are not considered. eg: scale.
    https://github.com/h5py/h5py/issues/428
        - strings must be fixed length for endianness interoperability
"""

from glob import glob
from math import sin, cos, radians
from os import remove, path
from struct import unpack

import numpy as np
import h5py as h5

from params import *
from params_base import *
from shared_bin import *

seis_file_list = glob(path.join(bin_output, '*_seis-?????.e3d'))

# python data-type format string
INT_S = 'i'
FLT_S = 'f'
# automatically swap bytes if required
if get_seis_swap(seis_file_list[0]):
    swapping_char = get_byteswap_char()
    INT_S = swapping_char + INT_S
    FLT_S = swapping_char + FLT_S

# read common data once
nt, dt, hh, rot = get_seis_common(seis_file_list[0], INT_S, FLT_S)
# rotation matrix for converting to 090, 000 (counter rot)
# ver is inverted (* -1) as E3D output is positive in down direction
theta = radians(rot)
rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
        [-sin(theta), -cos(theta), 0], \
        [0, 0, -1]])

# following need to be retrieved from params.py
mlat = float(MODEL_LAT)
mlon = float(MODEL_LON)
nx = int(nx)
ny = int(ny)
dx = int(dx_ts)
dy = int(dy_ts)
n_comps = N_MY_COMPS
comps = []
#for n, key in enumerate(MY_COMPS):
#    h5p.attrs['COMP_%d' % (n)] = np.string_(MY_COMPS[key])

# read every station in seis files
# keep virtual ones (name of 8 c-string bytes)
def process_seis_file(index):
    filename = seis_file_list[index]
    print('processing %s...' % (filename))
    fp = open(seis_file_list[index], 'rb')
    # speed optimise
    seek = fp.seek
    read = fp.read
    # shortcuts for processing input
    read_int = lambda : unpack(INT_S, read(SIZE_INT))[0]
    read_flt = lambda : unpack(FLT_S, read(SIZE_FLT))[0]

    # number of stations within this seis file
    num_stat = read_int()

    # map of array doesn't use RAM unless used
    comp_data = np.memmap(seis_file_list[index], \
            dtype = FLT_S, mode = 'r', \
            offset = SIZE_INT + num_stat * SIZE_SEISHEAD, \
            shape = (nt, num_stat, N_COMPS))

    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        stat = str(read(STAT_CHAR)).rstrip('\0')
        if len(stat) == VSTAT_LEN:
            # station is a virtual entry
            seek( - STAT_CHAR - 5 * SIZE_FLT - 4 * SIZE_INT, 1)
            x = read_int()
            y = read_int()

            # location in format (XXXXYYYY)
            g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))

            # grab station info from seis file
            #s_group.attrs['NAME'] = np.string_(stat)
            #s_group.attrs['X'] = x
            #s_group.attrs['Y'] = y
            #seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            #s_group.attrs['LAT'] = read_flt()
            #s_group.attrs['LON'] = read_flt()

            # read HF pairs
            try:
                hf0 = np.fromfile('%s_%s.%s' % \
                        (hf_prefix, g_name, MY_COMPS[0]), '>f4')
                hf1 = np.fromfile('%s_%s.%s' % \
                        (hf_prefix, g_name, MY_COMPS[1]), '>f4')
                hf2 = np.fromfile('%s_%s.%s' % \
                        (hf_prefix, g_name, MY_COMPS[2]), '>f4')
            except IOError:
                print('Warning: Cannot find HF files for %s.' % g_name)
                print('Skipping point %s.' % (g_name))
                continue

            # read LF pairs, rotate and reorientate
            lf = np.dot(np.dstack(( \
                    comp_data[..., stat_i, 0], \
                    comp_data[..., stat_i, 1], \
                    comp_data[..., stat_i, 2])), rot_matrix)[0]
            # convert to acceleration values
            lf0 = vel2acc(lf[..., 0])
            lf1 = vel2acc(lf[..., 1])
            lf2 = vel2acc(lf[..., 2])
            # apply amplification factors
            lf0 = ampdeamp(lf0, ampf, amp = True)
            lf1 = ampdeamp(lf1, ampf, amp = True)
            lf2 = ampdeamp(lf2, ampf, amp = True)
            # filter
            lf0 = bwfilter(lf0, 1.0, 'lowpass')
            lf1 = bwfilter(lf1, 1.0, 'lowpass')
            lf2 = bwfilter(lf2, 1.0, 'lowpass')

            # broadband join
            #
    fp.close()

# debug single threaded
for i in xrange(len(seis_file_list)):
    process_seis_file(i)

# multiprocessing
#p = Pool(procs)
#p.map(process_seis_file, seis_file_list)

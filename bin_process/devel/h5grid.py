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
from multiprocessing import Pool
from os import remove, path
from struct import unpack
import sys

import numpy as np
import h5py as h5

from params import *
from params_base import *
from shared_bin import *
from shared_ts import *

procs = 64
if len(sys.argv) > 1:
    procs = int(sys.argv[1])

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

# TODO: should have values for each location
vref = 500.0
vsite = 250.0
vpga = 500.0

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

    v_stat_info = []
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

            # read HF pairs
            try:
                hf = np.array([ \
                        np.fromfile('%s_%s.%s' % \
                                (hf_prefix, g_name, MY_COMPS[0]), '>f4'), \
                        np.fromfile('%s_%s.%s' % \
                                (hf_prefix, g_name, MY_COMPS[1]), '>f4'), \
                        np.fromfile('%s_%s.%s' % \
                                (hf_prefix, g_name, MY_COMPS[2]), '>f4')])
            except IOError:
                print('Warning: Cannot find HF files for %s.' % g_name)
                print('Skipping point %s.' % (g_name))
                continue

            pga = [0, 0, 0]
            ampf = [0, 0, 0]
            # process HF
            for tsi in xrange(len(hf)):
                # get pga
                pga[tsi] = np.max(np.abs(hf[tsi])) / 981.0
                # amplification factors
                ampf[tsi] = cb08_amp(dt, get_ft_len(nt), \
                        vref, vsite, vpga, pga[tsi])
                # apply amplification factors
                hf[tsi] = ampdeamp(hf[tsi], ampf, amp = True)
                # filter
                hf[tsi] = bwfilter(hf[tsi], dt, 1.0, 'highpass')

            # grab station info from seis file
            seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            v_stat_info.append({'NAME':stat, 'X':x, 'Y':y, \
                    'LAT':read_flt(), 'LON':read_flt(), \
                    'PGA_0':pga[0], 'PGA_1':pga[1], 'PGA_2':pga[2]})

            # read LF pairs, rotate and reorientate
            lf = np.dot(np.dstack(( \
                    comp_data[..., stat_i, 0], \
                    comp_data[..., stat_i, 1], \
                    comp_data[..., stat_i, 2])), rot_matrix)[0].T

            # difference in start timestep of LF, HF
            ddt = int(1.0 / dt)
            bb = np.empty(shape = (N_MY_COMPS, nt))
            # process LF, HF
            for tsi in xrange(len(lf)):
                # convert to acceleration
                lf[tsi] = vel2acc(lf[tsi], dt)
                # apply amplification factors
                lf[tsi] = ampdeamp(lf[tsi], ampf[tsi], amp = True)
                # filter
                lf[tsi] = bwfilter(lf[tsi], dt, 1.0, 'lowpass')

                # broadband
                bb[tsi] = np.cumsum((np.hstack(([0] * ddt, hf[tsi])) + \
                        np.hstack((lf[tsi], [0] * ddt))) * dt)
            # save broadband to file
            bb.T.astype(np.float32).tofile('bb_%s' % \
                    (g_name), format = '>f4')

    fp.close()
    return v_stat_info

# debug single threaded
seis_results = []
for i in xrange(len(seis_file_list)):
    seis_results.append(process_seis_file(i))

# multiprocessing
p = Pool(procs)
seis_results = p.map(process_seis_file, seis_file_list)

print('BroadBand Processing Complete')
print('Storing Results in HDF5')

h5p = h5.file('virtual.hdf5', 'w')
# add common metadata
#h5p.attrs['NSTAT'] = 
h5p.attrs['MLAT'] = mlat
h5p.attrs['MLON'] = mlon
h5p.attrs['NX'] = nx
h5p.attrs['NY'] = ny
h5p.attrs['DX'] = dx
h5p.attrs['DY'] = dy
h5p.attrs['NT'] = nt
h5p.attrs['DT'] = dt
h5p.attrs['HH'] = hh
h5p.attrs['ROT'] = rot
h5p.attrs['NCOMPS'] = N_MY_COMPS
for n, key in enumerate(MY_COMPS):
    h5p.attrs['COMP_%d' % (n)] = np.string_(MY_COMPS[key])

# add station data
for file_results in seis_results:
    for stat_i in file_results:
        g_name = '%s%s' % (str(stat_i['X']).zfill(4), str(stat_i['Y']).zfill(4))
        f_name = 'bb_%s' % (g_name)
        h5group = h5p.create_group(g_name)
        # have to save strings as fixed width because of h5py bug
        h5group.attrs['NAME'] = np.string_(stat_i['NAME'])
        for key in ['X', 'Y', 'LAT', 'LON', 'PGA_0', 'PGA_1', 'PGA_2']:
            h5group.attrs[key] = stat_i[key]
        # store data as components for each timestep
        h5group.create_dataset('VEL', (nt, N_MY_COMPS), \
                dtype = 'f', data = np.fromfile(f_name, dtype = 'f'))

h5p.close()

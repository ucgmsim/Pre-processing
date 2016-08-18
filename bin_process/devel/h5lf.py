#!/usr/bin/env python2
"""
Retrieves X, Y, Z components from SEIS files.
Saves to HDF5 file for all grid points.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 2 June 2016

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
from shared_bin import *

h5_file = 'virtual.hdf5'
h5p = h5.File(h5_file, 'w')
seis_file_list = glob(path.join(bin_output, '*_seis-?????.e3d'))

# python data-type format string
INT_S = 'i'
FLT_S = 'f'
# automatically swap bytes if required
if get_seis_swap(seis_file_list[0]):
    swapping_char = get_byteswap_char()
    INT_S = swapping_char + INT_S
    FLT_S = swapping_char + FLT_S

# read every station in seis files
# keep virtual ones (name of 8 c-string bytes)
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

    # number of stations within this seis file
    num_stat = read_int()

    # write common data to HDF5
    # only if this is the first seis file to be read
    try:
        h5p.attrs['NSTAT']
    # Fitzroy h5py throws TypeError (incorrectly)
    except (KeyError, TypeError):
        h5p.attrs['NSTAT'] = 0
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
        h5p.attrs['MLAT'] = float(MODEL_LAT)
        h5p.attrs['MLON'] = float(MODEL_LON)
        h5p.attrs['NX'] = int(nx)
        h5p.attrs['NY'] = int(ny)
        h5p.attrs['DX'] = int(dx_ts)
        h5p.attrs['DY'] = int(dy_ts)
        h5p.attrs['NT'] = nt
        h5p.attrs['DT'] = dt
        h5p.attrs['HH'] = hh
        h5p.attrs['ROT'] = rot
        h5p.attrs['NCOMPS'] = N_MY_COMPS
        for n, key in enumerate(MY_COMPS):
            h5p.attrs['COMP_%d' % (n)] = np.string_(MY_COMPS[key])


    # read at once, all component data below header as float array
    # major speedup compared to manual processing
    # TODO: look at reshaping comp_data to (NT, n_stat, N_BIN_COMPS)
    seek(SIZE_INT + num_stat * SIZE_SEISHEAD)
    comp_data = np.fromfile(fp, dtype = FLT_S)

    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        stat = str(read(STAT_CHAR)).rstrip('\0')
        if len(stat) == VSTAT_LEN:
            # station is a virtual entry
            h5p.attrs['NSTAT'] += 1

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

            # grab station info from seis file
            s_group.attrs['NAME'] = np.string_(stat)
            s_group.attrs['X'] = x
            s_group.attrs['Y'] = y
            seek(2 * SIZE_INT + 3 * SIZE_FLT, 1)
            s_group.attrs['LAT'] = read_flt()
            s_group.attrs['LON'] = read_flt()

            # prepare data for this station
            comps = np.dot(np.dstack(( \
                    comp_data[0 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[1 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[2 + stat_i * N_COMPS : : num_stat * N_COMPS])), rot_matrix)
            try:
                h5p.create_dataset(g_name + '/VEL', (nt, N_MY_COMPS), \
                        dtype='f', data = comps)
            except RuntimeError:
                print('failed to add station %s, replacing existing station' % (g_name))
                # this station had been added before?
                h5p[g_name + '/VEL'][...] = comps

    fp.close()

h5p.close()

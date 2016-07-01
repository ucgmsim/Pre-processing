#!/usr/bin/env python2
"""
Generates ASCII format version of seis file for each station/component.

@author Viktor Polak, Sung Bae
@date 5 April 2016

Replaces winbin-aio.csh. Re implemented in python. Using e3d.par.

USAGE: execute from current directory being the simulation directory containing the 'Vel' folder.
    if outside $sim_dir, link within: 'ln -s location/to/winbin-aio.py simdir/'

ISSUES: think of better filename

10 June 2016:
Include logic of fdbin2wcc and wcc_rotate binaries as python.
Removes dependencies of (very) slow and inefficient programs.
scale variable assumed to always be 1.0 (no reason to use other values?)
MODEL_ROT read directly from seis file, saves depending on variables
TSTRT assumed to be 0 (start at t=0), need to add functionality if needed
"""

from glob import glob
from math import sin, cos, radians
import os
import sys
sys.path.append(os.path.abspath(os.path.curdir))

import numpy as np

from shared import *
from shared_bin import *
from params import *

if len(sys.argv) > 1 and sys.argv[1] == 'test_mode':
    print('Running under test mode.')
    from postprocess_test.test_params import *

verify_files([FD_STATLIST])
verify_strings([run_name])
verify_user_dirs([vel_dir], reset = True) # prevent mixing with old files
verify_dirs([bin_output])

# stations to retrieve from seis files
stats = get_stations(FD_STATLIST)

filepattern = os.path.join(bin_output, run_name + '_seis*.e3d')
seis_file_list = glob(filepattern)

# python data-type format string
INT_S = 'i' 
FLT_S = 'f' 
if get_seis_swap(seis_file_list[0]):
    swapping_char = get_byteswap_char()
    INT_S = swapping_char + INT_S
    FLT_S = swapping_char + FLT_S

# read common data only from first file
fp = open(seis_file_list[0], 'rb')
read_flt = lambda : unpack(FLT_S, fp.read(SIZE_FLT))[0]
fp.seek(SIZE_INT * 5)
nt = unpack(INT_S, fp.read(SIZE_INT))[0]
dt = read_flt()
hh = read_flt()
rot = read_flt()
# rotation matrix for converting to 090, 000
# ver is inverted (* -1)
theta = radians(rot)
rot_matrix = np.array([[cos(theta), -sin(theta), 0], \
        [-sin(theta), -cos(theta), 0], \
        [0, 0, -1]])


for si, seis_file in enumerate(seis_file_list):
    print('processing input file %d of %d...' \
            % (si + 1, len(seis_file_list)))
    fp = open(seis_file, 'rb')
    fp_loaded = False
    # speed optimise
    seek = fp.seek
    read = fp.read

    num_stat = unpack(INT_S, fp.read(SIZE_INT))[0]

    for stat_i in xrange(num_stat):
        seek(SIZE_INT + (stat_i + 1) * SIZE_SEISHEAD - STAT_CHAR)
        stat = str(read(STAT_CHAR)).rstrip('\0')
        if stat in stats:
            stats.remove(stat)
            # station is of interest
            if not fp_loaded:
                # read at once, all component data below header as float array
                # major speedup compared to manual processing
                seek(SIZE_INT + num_stat * SIZE_SEISHEAD)
                comp_data = np.fromfile(fp, dtype = FLT_S)
                fp_loaded = True

            comps = np.dot(np.dstack(( \
                    comp_data[0 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[1 + stat_i * N_COMPS : : num_stat * N_COMPS], \
                    comp_data[2 + stat_i * N_COMPS : : num_stat * N_COMPS])), rot_matrix)[0]

            for i, comp in enumerate(MY_COMPS.values()):
                write_seis_ascii(vel_dir, comps[:, i], stat, comp, nt, dt)

    fp.close()

if len(stats) > 0:
    print('WARNING: the following stations were not found: %s' % str(stats))

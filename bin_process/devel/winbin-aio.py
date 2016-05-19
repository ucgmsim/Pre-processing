#!/usr/bin/env python2
"""
Generates ASCII format version of e3d file for each station/component.

@author Viktor Polak, Sung Bae
@date 5 April 2016

Replaces winbin-aio.csh. Re implemented in python. Using e3d.par.

USAGE: execute from current directory being the simulation directory containing the 'Vel' folder.
    if outside $sim_dir, link within: 'ln -s location/to/winbin-aio.py simdir/'

ISSUES: think of better filename
"""

import os
import sys
sys.path.append(os.path.abspath(os.path.curdir))
import shutil
from subprocess import call
from glob import glob

from shared import *
from params import *
if len(sys.argv) > 1 and sys.argv[1] == 'test_mode':
    print('Running under test mode.')
    from postprocess_test.test_params import *

fdbin2wcc_bin = os.path.join(wcc_prog_dir, 'fdbin2wcc')

verify_binaries([fdbin2wcc_bin])
verify_files([FD_STATLIST])
verify_logfiles([FILELIST])
verify_strings([output_prefix, scale, TSTRT, MODEL_ROT])
verify_user_dirs([vel_dir],reset=True) #Vel directory is best to start from empty
verify_user_dirs([bin_output])

STATS, LATS, LONS = get_stations(FD_STATLIST, True)
comp1 = str(int(90 + float(MODEL_ROT))).zfill(3)
comp2 = str(int(comp1) + 90).zfill(3)
COMPS = [comp1, comp2, 'ver']
AZ_COLUMN = 8

for stat in STATS:
    print(stat)


IX = ['0', '1', '2']
FLIP = ['1', '1', '-1']

list_handle = open(FILELIST, 'w')

filepattern = os.path.join(bin_output, output_prefix+ '_seis*.e3d')
print filepattern
list_of_files = '\n'.join([file_path for file_path in glob(filepattern)]) + '\n'
print list_of_files
list_handle.write(list_of_files)
list_handle.close()

for s_index, stat in enumerate(STATS):
    print(LONS[s_index] + ' ' + LATS[s_index] + ' ' + stat)

    statfile = os.path.join(vel_dir, stat)

    for c_index, comp in enumerate(COMPS):
        cmd = [fdbin2wcc_bin, 'all_in_one=1', 'filelist=' + FILELIST,
                'ix=' + IX[c_index], 'scale=' + scale, 'flip=' + FLIP[c_index], 'tst=' + TSTRT,
                'stat=' + stat, 'comp=' + comp, 'outfile=' + statfile + '.' + comp,
                'nseis_comps=9', 'swap_bytes=0', 'tst2zero=0', 'outbin=0'] 
        print ' '.join(cmd)
        call(cmd)

    # XXX: THIS IS WRONG. ROT IS RELATED TO MODEL_ROT. FORWARDS OR BACKWARDS?
    # ONLY CALL THIS IF MODEL_ROT != 0?
    cmd = [os.path.join(wcc_prog_dir,'wcc_rotate'), \
        'filein1=' + statfile + '.' + COMPS[0],'inbin1=0', \
        'filein2=' +  statfile + '.' + COMPS[1],'inbin2=0', \
        'fileout1=' + statfile + '.000', 'outbin1=0', \
        'fileout2=' + statfile + '.090', 'outbin2=0', \
        'rot=0.0']
    print ' '.join(cmd)
    call(cmd)

    # only getting rid of pre-rotated files
    if int(float(MODEL_ROT)) != 0:
        for i in xrange(len(COMPS) - 1):
            try:
                os.remove('%s.%s'%(statfile,COMPS[i]))
            except OSError:
                pass



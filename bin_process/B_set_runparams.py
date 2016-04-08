#!/usr/bin/env python2

from __future__ import print_function
import fileinput
import sys
from shutil import copyfile

from A_edit_runparams import *

# note this only works with setting strings (value will be surrounded by '')
def change_var(filename, variable, value):
    for line in fileinput.input(filename, inplace = True):
        if line.startswith(variable + ' = '):
            line = variable + ' = \'' + value + '\'\n'
        # don't include extra '\n'
        sys.stdout.write(line)

try:
    copyfile(DEFAULT_PARFILE, PARFILE)
except IOError:
    print('Cannot copy DEFAULT_PARFILE to PARFILE! in ' + __file__, file = sys.stderr)
    raise

try:
    # default parfile is appended to
    par_handle = open(PARFILE, 'a')
except IOError:
    print('PARFILE cannot be opened to append data in ' + __file__, file = sys.stderr)
    raise

configs = ['version=' + VERSION + '-mpi', \
'name=' + RUN_NAME, \
'nx=' + NX, \
'ny=' + NY, \
'nz=' + NZ, \
'h=' + HH, \
'nt=' + NT, \
'dt=' + DT, \
 \
'bfilt=4', \
'flo=' + FLO, \
'fhi=0.0', \
'bforce=0', \
'pointmt=0', \
'dblcpl=0', \
'ffault=2', \
'faultfile=' + SRF_FILE, \
 \
'model_style=1', \
# only for the 1D velocity model
#'model=' + FD_VMODFILE, \
'vmoddir=' + VMODDIR, \
'pmodfile=' + PMOD, \
'smodfile=' + SMOD, \
'dmodfile=' + DMOD, \
'qpfrac=100', \
'qsfrac=50', \
'qpqs_factor=2.0', \
'fmax=25.0', \
'fmin=0.01', \
'vmodel_swapb=1', \
 \
'modellon=' + MODEL_LON, \
'modellat=' + MODEL_LAT, \
'modelrot=' + MODEL_ROT, \
 \
'enable_output_dump=1', \
'dump_itinc=' + DUMP_ITINC, \
'main_dump_dir=' + MAIN_OUTPDIR, \
'nseis=1', \
'seiscords=' + STATCORDS, \
'seisdir=' + TMP_SEISDIR, \
 \
'ts_xy=1', \
'iz_ts=1', \
'dtts=' + DT_TS, \
'dxts=' + DX_TS, \
'dyts=' + DY_TS, \
'dzts=' + DZ_TS, \
 \
'enable_restart=' + ENABLE_RESTART, \
'restartdir=' + MAIN_RESTARTDIR, \
'restart_itinc=' + RESTART_ITINC, \
'read_restart=' + READ_RESTART, \
'restartname=' + RUN_NAME]

par_handle.write('\n'.join(configs) + '\n')
par_handle.close()

############# updating variables in other files #############
change_var('winbin-aio.py', 'RUN', RUN_NAME)
change_var('merge_tsP3.py', 'FILEROOT', RUN_NAME)
change_var('gen_ts_default.py', 'FILEROOT', RUN_NAME)
change_var('plot_ts_bluefern.py', 'NAME', RUN_NAME)


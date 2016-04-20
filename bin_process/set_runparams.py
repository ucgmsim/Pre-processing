#!/usr/bin/env python2

from __future__ import print_function
import fileinput
import sys
from shutil import copyfile
from params import *


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

configs = ['version=' + version + '-mpi', \
'name=' + run_name, \
'extended_name=' + extended_run_name, \
'nproc=' + n_proc, \
'nx=' + nx, \
'ny=' + ny, \
'nz=' + nz, \
'h=' + hh, \
'nt=' + nt, \
'dt=' + dt, \
 \
'bfilt=4', \
'flo=' + flo, \
'fhi=0.0', \
'bforce=0', \
'pointmt=0', \
'dblcpl=0', \
'ffault=2', \
'faultfile=' + srf_file, \
 \
'model_style=1', \
# only for the 1D velocity model
#'model=' + FD_VMODFILE, \
'vmoddir=' + vel_mod_dir, \
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
'main_dump_dir=' + bin_output, \
'nseis=1', \
'seiscords=' + stat_coords, \
'seisdir=' + seis_tmp_dir, \
 \
'ts_xy=1', \
'iz_ts=1', \
'dtts=' + dt_ts, \
'dxts=' + dx_ts, \
'dyts=' + dy_ts, \
'dzts=' + dz_ts, \
'ts_start=' + ts_start, \
'ts_inc=' + ts_inc, \
'ts_total=' + ts_total, \
'swap_bytes=' + swap_bytes, \
'lonlat_out=' + lonlat_out, \
'scale=' + scale, \
 \
'enable_restart=' + ENABLE_RESTART, \
'restartdir=' + restart_dir, \
'restart_itinc=' + RESTART_ITINC, \
'read_restart=' + READ_RESTART, \
'restartname=' + RUN_NAME, \
\
# plot_ts.sh
'plot_main_title="' + plot_main_title + '"', \
'plot_sub_title="' + plot_sub_title + '"', \
'plot_option=' + plot_option, \
\
# other locations
'wcc_prog_dir="' + wcc_prog_dir + '"', \
'vel_mod_params_dir="' + vel_mod_params_dir + '"', \
'stat_file="' + stat_file + '"']

par_handle.write('\n'.join(configs) + '\n')
par_handle.close()



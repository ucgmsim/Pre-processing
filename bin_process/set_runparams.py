#!/usr/bin/env python2

from __future__ import print_function
import fileinput
import sys
import os.path
sys.path.append(os.path.abspath(os.path.curdir))
from shutil import copyfile
from params import *

try:
    copyfile(default_parfile, parfile)
except IOError:
    print('Cannot copy default_parfile to parfile! in ' + __file__, file = sys.stderr)
    raise

try:
    # default parfile is appended to
    par_handle = open(parfile, 'a')
except IOError:
    print('parfile cannot be opened to append data in ' + __file__, file = sys.stderr)
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
'ts_file=' + ts_file, \
'ts_out_dir="' + ts_out_dir + '"', \
'ts_out_prefix="' + ts_out_prefix + '"', \
'swap_bytes=' + swap_bytes, \
'lonlat_out=' + lonlat_out, \
'scale=' + scale, \
 \
'enable_restart=' + ENABLE_RESTART, \
'restartdir=' + restart_dir, \
'restart_itinc=' + RESTART_ITINC, \
'read_restart=' + READ_RESTART, \
'restartname=' + run_name, \
\
# plot_ts.sh
'plot_main_title="' + plot_main_title + '"', \
'plot_sub_title="' + plot_sub_title + '"', \
'plot_option=' + plot_option, \
'plot_palette=' + plot_palette, \
'plot_sites=' + plot_sites, \
'plot_s_pos=' + plot_s_pos, \
'plot_s_lon=' + plot_s_lon, \
'plot_s_lat=' + plot_s_lat, \
'plot_s_sym="' + plot_s_sym + '"', \
'plot_s_fil="' + plot_s_fil + '"', \
'plot_s_lin="' + plot_s_lin + '"', \
'plot_x_org=' + plot_x_org, \
'plot_y_org=' + plot_y_org, \
'plot_x_inch=' + plot_x_inch, \
'plot_x_shift=' + plot_x_shift, \
'plot_x_min=' + plot_x_min, \
'plot_x_max=' + plot_x_max, \
'plot_y_min=' + plot_y_min, \
'plot_y_max=' + plot_y_max, \
'plot_region=' + plot_region, \
'plot_ts_x_min=' + plot_ts_x_min, \
'plot_ts_x_max=' + plot_ts_x_max, \
'plot_ts_y_min=' + plot_ts_y_min, \
'plot_ts_y_max=' + plot_ts_y_max, \
'plot_ts_region=' + plot_ts_region, \
'plot_dx=' + plot_dx, \
'plot_dy=' + plot_dy, \
'plot_ps_dir="' + plot_ps_dir + '"', \
'plot_png_dir="' + plot_png_dir + '"', \
'plot_res=' + plot_res, \
'plot_orig_dt=' + plot_orig_dt, \
'plot_comps=' + plot_comps, \
'plot_topo_file="' + plot_topo_file + '"', \
'plot_topo_illu="' + plot_topo_illu + '"', \
'plot_topo_a_min=' + plot_topo_a_min, \
'plot_topo_a_inc=' + plot_topo_a_inc, \
'plot_topo_a_max=' + plot_topo_a_max, \
'plot_topo_a_below=' + plot_topo_a_below, \
'plot_fault_add_plane="' + plot_fault_add_plane + '"', \
'plot_fault_line=' + plot_fault_line, \
'plot_fault_top_edge=' + plot_fault_top_edge, \
'plot_fault_hyp_open=' + plot_fault_hyp_open, \
\
# other locations
'wcc_prog_dir="' + wcc_prog_dir + '"', \
'vel_mod_params_dir="' + vel_mod_params_dir + '"', \
'global_root="' + global_root + '"', \
'sim_dir="' + sim_dir + '"', \
'fault_file="' + fault_file + '"', \
'stat_file="' + stat_file + '"']

par_handle.write('\n'.join(configs) + '\n')
par_handle.close()



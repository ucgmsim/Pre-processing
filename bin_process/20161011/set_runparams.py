#!/usr/bin/env python2
"""
Generates 'e3d.par' from the default set, appending new key value pairs of parameters.

@author Viktor Polak, Sung Bae
@date 6 April 2016

Replaces set_runparams.csh. Converted to python, handles params set in separate file.

USAGE: edit params.py (not this file), execute this file from the same directory.

ISSUES: remove default values in e3d_default.par where not needed.
"""

from __future__ import print_function
import fileinput
import sys
import os.path
sys.path.append(os.path.abspath(os.path.curdir))
from shutil import copyfile
import shared

params_uncertain = 'params_uncertain.py'


if __name__ == '__main__' :

    # attempt to append template file before importing params
    try:
        # throws NameError if var not set, AssertionError if blank
        assert(params_override != '')
        # copy to temp file
        copyfile('params.py', 'params_joined.py')
        # append to temp
        with open('params_joined.py', 'a') as fp:
            with open('params_override_' + params_override + '.py', 'r') as tp:
                fp.write(tp.readlines())
        # import temp
        import params_joined
        os.remove('params_joined.py')
    except (AssertionError, NameError, ImportError, OSError):
        from params import *


    p1={}
    p2=[]  
    for i, srf_file in enumerate(srf_files):
        srf_file_basename = os.path.basename(srf_file).split('.')[0] #take the filename only
        p1['lf_sim_dir'] = os.path.join(lf_sim_root_dir,srf_file_basename)
        shared.verify_user_dirs([p1['lf_sim_dir']])

        p1['restart_dir'] = os.path.join(p1['lf_sim_dir'], 'Restart')
        p1['bin_output'] = os.path.join(p1['lf_sim_dir'], 'OutBin')
        p1['SEISDIR'] = p1['bin_output']
        p1['ts_file'] = os.path.join(p1['bin_output'], run_name+ '_xyts.e3d') #the file created by merge_ts

        p1['log_dir'] = os.path.join(p1['lf_sim_dir'], 'Rlog')
        p1['slipout_dir'] = os.path.join(p1['lf_sim_dir'],'SlipOut')
        p1['vel_dir'] = os.path.join(p1['lf_sim_dir'],'Vel')
        p1['t_slice_dir'] = os.path.join(p1['lf_sim_dir'], 'TSlice')
         # output dirs and resolution (dpi)
        p1['plot_ps_dir'] = os.path.join(p1['t_slice_dir'], 'PlotFiles') #only written to e3d.par
        p1['plot_png_dir'] = os.path.join(p1['t_slice_dir'], 'Png') #only written to e3d.par

    
        p1['ts_out_dir'] = os.path.join(p1['t_slice_dir'], 'TSFiles')
        p1['ts_out_prefix'] = os.path.join(p1['ts_out_dir'], run_name)
        p1['FILELIST'] = os.path.join(p1['lf_sim_dir'],'fdb.filelist')

        write_to_py(os.path.join(p1['lf_sim_dir'],params_uncertain),p1)

        
        p2['version']= version+'-mpi'
        p2['name']= run_name
        p2['nproc']=n_proc

        p2['nx']= nx
        p2['ny']= ny
        p2['nz']= nz
        p2['h']= hh
        p2['nt']= nt
        p2['dt']= dt
        p2['bfilt']=4
        p2['flo'] = flo
        p2['fhi']=0.0
        
        p2['bforce']=0
        p2['pointmt']=0
        p2['dblcpl']=0
        p2['ffault']=2
        p2['faultfile']=srf_file

        p2['model_style']=1
        # only for the 1D velocity model
        #'model=' + FD_VMODFILE, \
        p2['vmoddir']= vel_mod_dir
        p2['pmodfile']= PMOD
        p2['smodfile']= SMOD
        p2['dmodfile']= DMOD
        p2['qpfrac']=100
        p2['qsfrac']=50
        p2['qpqs_factor']=2.0
        p2['fmax']=25.0
        p2['fmin']=0.01
        p2['vmodel_swapb']=1

        p2['modellon']= MODEL_LON
        p2['modellat']= MODEL_LAT
        p2['modelrot']= MODEL_ROT
        
        p2['enable_output_dump']=1
        p2['dump_itinc']=DUMP_ITINC
        p2['main_dump_dir']= p1['bin_output']
        p2['nseis']=1
        p2['seiscords']= stat_coords
        p2['seisdir']= seis_tmp_dir
        p2['ts_xy']= ts_xy
        p2['iz_ts']= iz_ts
        p2['ts_xz']= ts_xz
        p2['iy_ts']= iy_ts
        p2['ts_yz']= ts_yz
        p2['ix_ts']= ix_ts
        p2['dtts']= dt_ts
        p2['dxts']= dx_ts
        p2['dyts']= dy_ts
        p2['dzts']= dz_ts
        p2['ts_start']= ts_start
        p2['ts_inc']= ts_inc
        p2['ts_total']= ts_total
        p2['ts_file']= p1['ts_file']
        p2['ts_out_dir']= p1['ts_out_dir']
        p2['ts_out_prefix']= p1['ts_out_prefix']
        p2['swap_bytes']= swap_bytes
        p2['lonlat_out']= lonlat_out
        p2['scale']= scale
        p2['enable_restart']= ENABLE_RESTART
        p2['restartdir']= p1['restart_dir']
        p2['restart_itinc']= RESTART_ITINC
        p2['read_restart']= READ_RESTART
        p2['restartname']= run_name
        p2['logdir']= p1['log_dir']
        p2['slipout']= p1['slipout_dir']+'/slipout-k2'
        # extras found in default parfile
        p2['span']=1
        p2['intmem']=1
        p2['maxmem']=1500
        p2['order']=4
        p2['model_style']=1
        p2['elas_only']=0
        p2['freesurf']=1
        p2['dampwidth']=0
        p2['qbndmax']=100.0
        p2['stype']='2tri-p10-h20'
        p2['tzero']=0.6
        p2['geoproj']=1
        p2['report']=100
        p2['all_in_one']=1
        p2['xseis']=0
        p2['iy_xs']=60
        p2['iz_xs']=1
        p2['dxout']=1
        p2['yseis']=0
        p2['ix_ys']=100
        p2['iz_ys']=1
        p2['dyout']=1
        p2['zseis']=0
        p2['ix_zs']=100
        p2['iy_zs']=50
        p2['dzout']=1
        p2['ts_xz']=0
        p2['iy_ts']=2
        p2['ts_yz']=0
        p2['ix_ts']=99
        # plot_ts.sh
        p2['plot_main_title']= plot_main_title
        p2['plot_sub_title']= plot_sub_title
        p2['plot_option']= plot_option
        p2['plot_palette']= plot_palette
        p2['plot_sites']= plot_sites
        p2['plot_s_pos']= plot_s_pos
        p2['plot_s_lon']= plot_s_lon
        p2['plot_s_lat']= plot_s_lat
        p2['plot_s_sym']= plot_s_sym
        p2['plot_s_fil']= plot_s_fil 
        p2['plot_s_lin']= plot_s_lin
        p2['plot_x_org']= plot_x_org
        p2['plot_y_org']= plot_y_org
        p2['plot_x_inch']= plot_x_inch
        p2['plot_x_shift']= plot_x_shift
        p2['plot_x_min']= plot_x_min
        p2['plot_x_max']= plot_x_max
        p2['plot_y_min']= plot_y_min
        p2['plot_y_max']= plot_y_max
        p2['plot_region']= plot_region
        p2['plot_ts_x_min']= plot_ts_x_min
        p2['plot_ts_x_max']= plot_ts_x_max
        p2['plot_ts_y_min']= plot_ts_y_min
        p2['plot_ts_y_max']= plot_ts_y_max
        p2['plot_ts_region']= plot_ts_region
        p2['plot_dx']= plot_dx
        p2['plot_dy']= plot_dy
        p2['plot_ps_dir']= p1['plot_ps_dir']
        p2['plot_png_dir']= p1['plot_png_dir']
        p2['plot_res']= plot_res
        p2['plot_orig_dt']= plot_orig_dt
        p2['plot_comps']= plot_comps
        p2['plot_topo_file']= plot_topo_file
        p2['plot_topo_illu']= plot_topo_illu 
        p2['plot_topo_a_min']= plot_topo_a_min
        p2['plot_topo_a_inc']= plot_topo_a_inc
        p2['plot_topo_a_max']= plot_topo_a_max
        p2['plot_topo_a_below']= plot_topo_a_below
        p2['plot_fault_add_plane']= plot_fault_add_plane
        p2['plot_fault_line']= plot_fault_line
        p2['plot_fault_top_edge']= plot_fault_top_edge
        p2['plot_fault_hyp_open']= plot_fault_hyp_open
        # other locations
        p2['wcc_prog_dir']= wcc_prog_dir
        p2['vel_mod_params_dir']= vel_mod_params_dir
        p2['global_root']= global_root
        p2['sim_dir']= sim_dir 
        p2['fault_file']= fault_file 
        p2['stat_file']= stat_file 
        p2['grid_file']= GRIDFILE
        p2['model_params']= MODELPARAMS
        
        write_to_py(os.path.join(p1['lf_sim_dir'],parfile),p)


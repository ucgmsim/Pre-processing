import os.path
if __name__ == '__main__':
    # in case someone has passed this file as a parameter to python
    print('\nDo not run ' + __file__ + \
            '\nOnly edit variables in ' + __file__ + ' and then run the set param script.\n')
    exit()


# parameters are written to files, need to be strings

######## Global CONSTANTS ########

run_name = '2011Feb22_m6pt2bev01_Cantv1.64'
version = '3.0.4'
extended_run_name = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_Test'

# keep as int, processed before writing to file
n_proc_x = 16
n_proc_y = 16
n_proc_z = 8


flo = '1.0'
# spacial grid spacing (km)
hh = '0.100'
# x, y, z distance (km)
nx = '1400'
ny = '1200'
nz = '460'
# number of timesteps n, time step Dt
nt = '20000'
dt = '0.005'

MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'

DUMP_ITINC = '4000'

dt_ts = '20'
# x, y, z decimation along axis
dx_ts = '5'
dy_ts = '5'
# 1 for xy time slice
dz_ts = '1'

# which time slices to iterate over
ts_start = '0'     # first one is 0
ts_inc = '1'       # increment, larger than 1 to skip
ts_total = '400'   # number of slices to generate. sim time = ts_total * ORIG_DT

# swap_bytes 0/1 no/yes - should be 1 if
#   ts_file created on supercomp and this file is run on laptop; zero if run within supercomputer)
# the three lines below not used if the TSFiles are created using 'gen_ts.py' on supercomputer
#   and copied to local computer beforehand.
swap_bytes = '0'
lonlat_out = '1'
scale = '1.0'

# TODO: remove default? would simplify, remove duplicates. anything lost (is appended to)?
default_parfile = 'e3d_default.par'
parfile = 'e3d.par'

#FD_VMODFILE = 'Cant1D_v1.fd_modfile'     #This line was for a 1D Vp,Vs,rho model
# set names of P,S,D files
PMOD = 'vp3dfile.p'
SMOD = 'vs3dfile.s'
DMOD = 'rho3dfile.d'

v_mod_ver = 'v1.64'

ENABLE_RESTART = '0'
READ_RESTART = '0'
RESTART_ITINC = '20000'


######## Variables derived from Global CONSTANTS ##########

n_proc = str(n_proc_x * n_proc_y * n_proc_z)

# works on Windows and POSIX paths
system_root = os.path.abspath(os.sep)
user_scratch = os.path.join(system_root, 'hpc', 'scratch', os.getenv('USER'))

# things that everyone doesn't have, eg. binaries are within here
global_root = os.path.join(system_root, 'hpc', 'scratch', 'nesi00213')
# things that people have their own copies of, eg. RunFolder
user_root = os.path.expanduser('~')

# directories - main. change global_root with user_root as required
run_dir = os.path.join(user_root, 'RunFolder')
srf_dir = os.path.join(user_root, 'RupModel')
stat_dir = os.path.join(user_root, 'StationInfo')
vel_mod_params_dir = os.path.join(user_root, 'VelocityModel', 'ModelParams')
vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel', v_mod_ver)
wcc_prog_dir = os.path.join(global_root, 'EMOD3D', 'WccFormat', 'bin', 'ppc64-Linux-p2n14-c-gcc')
seis_tmp_dir = os.path.join(user_scratch, extended_run_name, 'SeismoBin')
# directories - derived
sim_dir = os.path.join(run_dir, extended_run_name)
restart_dir = os.path.join(sim_dir, 'Restart')
bin_output = os.path.join(sim_dir, 'OutBin')

# files
srf_file = os.path.join(srf_dir, 'Srf', 'm6.20-16.0x9.0_s560.srf')
stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_nz01-h0.100.statcords')



############# winbin-aio ##############
# -------- START OF USER DEFINED INPUT -------------


# Define location of FD output files to use
SEISDIR = bin_output

# Define folder to write velocity seismograms to
VELDIR = os.path.join(sim_dir,'Vel')

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'
SCALE = '1.0'

# Define the directory with the station information, and components
FD_STATLIST = stat_dir + '/fd_nz01-h0.100.ll'


############### gen_ts ###################

ts_file = os.path.join(bin_output, extended_run_name + '_xyts.e3d') #the file created by merge_ts

# TODO: not used anywhere???
# ABSMAX =1 vector magnitude of all 3 components is output into <fileroot>.0
#        =0 3 components are output respectively into <fileroot>.0, <fileroot>.1, <fileroot>.2
ABSMAX = '1'

LONLAT_OUT = '1' # LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice
GRIDFILE = vel_mod_params_dir + '/gridout_nz01-h' + hh # GRIDFILE is the file containing the local (x,y,z) coordinates for this 3D run

t_slice_dir = os.path.join(sim_dir, 'TSlice')
ts_out_dir = os.path.join(t_slice_dir, 'TSFiles')
ts_out_prefix = os.path.join(ts_out_dir, run_name)

################## plot_ts ####################


# main title and subtitles that appears at the top of the picture
plot_main_title = 'Mw6.2 22 Feb 2011 Earthquake'
plot_sub_title = 'Beavan 1Fault, Stoch Slip, Chch 1D VM'

# input can be provided either via the "._xyts.e3d" file (option=1)
# or via all of the individual TSFiles, which are obtained from gen_ts.csh 
# within BlueFern and then scp'd to the local machine (option=2)
plot_option = '2'

# information required for option 2.
# ----------------------------------

# details of the spatial and termporal discretization and spacing
plot_orig_dt = '0.1'  # time step of the time slice output (this is DT*DT_TS from the run files)


# swap_bytes 0/1 no/yes - should be 1 if
#   ts_file created on supercomp and this file is run on laptop; zero if run within supercomputer)
# the three lines below not used if the TSFiles are created using 'gen_ts.py' on supercomputer
#   and copied to local computer beforehand.
swap_bytes = '0'
lonlat_out = '1'
scale = '1.0'

#components to consider (0= when only looking at vector magnitude [i.e. ABSMAX=1 further down].
plot_comps = '( 0 )'


# output dirs and resolution (dpi)
plot_ps_dir = os.path.join(t_slice_dir, 'PlotFiles')
plot_png_dir = os.path.join(t_slice_dir, 'Png')
plot_res = '140'   # 720 for PDF quality

TOPODIR = global_root + '/PlottingData/TopoData'
TOPO_FILE = TOPODIR + '/srtm_71_21.grd'
ILLU = '-I' + TOPODIR + '/srtm_71_21_i5.grd'
TOPO_FILE2 = TOPODIR + '/etopo2.grd'

# set PALETTE = '-Crelief.cpt'
plot_palette = '-Cgray'

# specific locations to display on figure (not currently used)
plot_sites = '(Rolleston Darfield Lyttelton Akaroa Kaiapoi Rakaia Oxford)'
# alignment; 2letter, L,C,R (left, center, right); T,M,B (top, middle, bottom)
plot_s_pos = "(RB CB LM RB LB RT LB)"
plot_s_lon = "(172.3791667 172.1116667 172.7194444 172.9683333 \
         172.6569444 172.0230556 172.1938889)"
plot_s_lat = "(-43.59083333 -43.48972222 -43.60305556 -43.80361111 \
         -43.38277778 -43.75611111 -43.29555556)"
# specifying plotting preferences for site locations
plot_s_sym = 'c0.10'         # symbol
plot_s_fil = '220/220/220'   # fill color
plot_s_lin = '1,000/000/000' # line?

# location of offset plotting (for when 3 component plotting used - not currently utilized)
plot_x_org = '1.15'
plot_y_org = '2.5'
plot_x_inch = '5.0'
plot_x_shift = '%0.6f' % (float(plot_x_inch) + 0.2)

# specify the maximum plotting window for the basemap
plot_x_min = '171.75'
plot_x_max = '173.00'
plot_y_min = '-44.00'
plot_y_max = '-43.20'
# region in GMT format
plot_region = '/'.join([plot_x_min, plot_x_max, plot_y_min, plot_y_max])

# specify the plotting region for the time slice

plot_ts_x_min = '171.75'
plot_ts_x_max = '173.00'
plot_ts_y_min = '-44.00'
plot_ts_y_max = '-43.20'
# region in GMT format
plot_ts_region = '/'.join([plot_ts_x_min, plot_ts_x_max, plot_ts_y_min, plot_ts_y_max])
# specify the increments of X/Y (cartesian coords) for masks etc.
plot_dx = '0.002'
plot_dy = '0.002'


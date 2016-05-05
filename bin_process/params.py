import os.path
if __name__ == '__main__':
    # in case someone has passed this file as a parameter to python
    print('\nDo not run ' + __file__ + \
            '\nOnly edit variables in ' + __file__ + ' and then run the set param script.\n')
    exit()

params_override = '2010Sep4'

# parameters are written to files, need to be strings

######## Global CONSTANTS ########
run_name = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_Test'
version = '3.0.4'
# directory name inder RunFolder
extended_run_name = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_Test'
# prefix used for naming OutBin/prefix{_xyts.e3d,_seis.e3d}
output_prefix = run_name

# keep as int, processed before writing to file
# only product is stored
n_proc_x = 8
n_proc_y = 8
n_proc_z = 8


flo = '1.0'
# spacial grid spacing (km)
hh = '0.100'
# x, y, z distance (km)
nx = '1400'
ny = '1200'
nz = '460'
# cap number of timesteps in simulation, not all timesteps have outputs
# max simulation timeperiod = nt * dt eg: 10,000 * 0.005 = 50 seconds
nt = '20000'
# dt should be 0.005 (or smaller), small increments required in simulation
dt = '0.005'

MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'

# how often to save outputs (measured in simulation timesteps)
DUMP_ITINC = '4000'

# output timestep in multiples of simulation timestep
# eg: simulation dt 0.005 sec * dt_ts 20 = 0.1 second increments
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


n_proc = str(n_proc_x * n_proc_y * n_proc_z)

# things that everyone doesn't have, eg. binaries are within here
global_root = '/nesi/projects/nesi00213'


# things that people have their own copies of, eg. RunFolder
user_root = global_root   # os.path.expanduser('~')   #NOT recommended to use separate user_root in Fitzroy.

# works on Windows and POSIX paths
user_scratch = os.path.join(user_root, 'scratch', os.getenv('USER'))

# directories - main. change global_root with user_root as required
run_dir = os.path.join(user_root, 'RunFolder')
srf_dir = os.path.join(user_root, 'RupModel')
stat_dir = os.path.join(user_root, 'StationInfo')
vel_mod_params_dir = os.path.join(user_root, 'VelocityModel/ModelParams')

vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel', v_mod_ver)
wcc_prog_dir = os.path.join(global_root, 'EMOD3D/WccFormat/bin/powerpc-AIX-nesi2-xlc')

seis_tmp_dir = os.path.join(user_scratch, extended_run_name, 'SeismoBin')

# directories - derived
sim_dir = os.path.join(run_dir, extended_run_name)
restart_dir = os.path.join(sim_dir, 'Restart')
bin_output = os.path.join(sim_dir, 'OutBin')

# files
srf_file = os.path.join(srf_dir, '2010Sept4_m7pt1/Srf/bev01.srf')
stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_nz01-h0.100.statcords')



############# winbin-aio ##############
# -------- START OF USER DEFINED INPUT -------------


# Define location of FD output files to use
SEISDIR = bin_output

# Define folder to write velocity seismograms to
vel_dir = os.path.join(sim_dir,'Vel')

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'
SCALE = '1.0'

# Define the directory with the station information, and components
FD_STATLIST = stat_dir + '/fd_nz01-h0.100.ll'
print FD_STATLIST

############### gen_ts ###################

#SIMDIR = MAIN_OUTPDIR #the directory name where the 3D FD run output exists

ts_file = os.path.join(bin_output, output_prefix+ '_xyts.e3d') #the file created by merge_ts
HH = '0.100' # HH is the grid spacing used for the 3D FD run
DXTS = '5' # DXTS is the decimation along the X-axis for the time slice output (same as dxts in e3d.par)
DYTS = '5' # DYTS is the decimation along the Y-axis for the time slice output
DZTS = '1' # DZTS is the decimation along the Z-axis for the time slice output (=1 for xy time slice)

#FD_STATLIST = stat_dir + '/fd_nz01-h0.100.ll'

# TODO: not used anywhere???
# ABSMAX =1 vector magnitude of all 3 components is output into <fileroot>.0
#        =0 3 components are output respectively into <fileroot>.0, <fileroot>.1, <fileroot>.2
ABSMAX = '1'

TS_INC = 1 # TS_INC is the increment for output time slices (=2 means every other slice is output)
TS_START = 0 # TS_START is the starting time slice (0 is first one)
TS_TOTAL = 400 # the total number of slices to generate. the total sim time = TS_TOTAL * ORIG_DT

LONLAT_OUT = '1' # LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice

GRIDFILE = vel_mod_params_dir + '/gridout_nz01-h' + hh # GRIDFILE is the file containing the local (x,y,z) coordinates for this 3D run

t_slice_dir = os.path.join(sim_dir, 'TSlice')
ts_out_dir = os.path.join(t_slice_dir, 'TSFiles')
ts_out_prefix = os.path.join(ts_out_dir, run_name)

# swap_bytes 0/1 no/yes - should be 1 if
#   ts_file created on supercomp and this file is run on laptop; zero if run within supercomputer)
# the three lines below not used if the TSFiles are created using 'gen_ts.py' on supercomputer
#   and copied to local computer beforehand.
swap_bytes = '0'
lonlat_out = '1'
scale = '1.0'

################## plot_ts ####################


# main title and subtitles that appears at the top of the picture
plot_main_title = 'Mw7.1 4 Sept 2010 Earthquake'
plot_sub_title = 'Beavan 1Fault, Stoch Slip, v1.64'

# input can be provided either via the "._xyts.e3d" file (option=1)
# or via all of the individual TSFiles, which are obtained from gen_ts.csh 
# within BlueFern and then scp'd to the local machine (option=2)
plot_option = '2'

# option 2 settings follow
# ----------------------------------

# details of the spatial and termporal discretization and spacing
plot_orig_dt = '0.1'  # time step of the time slice output (this is DT*DT_TS from the run files)



#components to consider (0= when only looking at vector magnitude [i.e. ABSMAX=1 further down].

plot_comps = '( 0 )'


# output dirs and resolution (dpi)
plot_ps_dir = os.path.join(t_slice_dir, 'PlotFiles')
plot_png_dir = os.path.join(t_slice_dir, 'Png')
plot_res = '140'   # 720 for PDF quality

# topography
plot_topo_dir = os.path.join(global_root, 'PlottingData/TopoData')
plot_topo_file = os.path.join(plot_topo_dir, 'srtm_71_21.grd')
plot_topo_illu = '-I' + os.path.join(plot_topo_dir, 'srtm_71_21_i5.grd')
plot_topo_file_2 = os.path.join(plot_topo_dir, 'etopo2.grd')
# altitude in killometers? for topo
plot_topo_a_min = '2'
plot_topo_a_inc = '10'
plot_topo_a_max = '80'
plot_topo_a_below = 'NaN'

# fault plane resources
fault_plane_dir = os.path.join(global_root, 'PlottingData', 'sourcesAndStrongMotionStations')
fault_file = os.path.join(fault_plane_dir, 'bev01_DarfieldFaultPlane.xy')
plot_fault_add_plane = os.path.join(fault_plane_dir, 'addStandardFaultPlane.sh')
plot_fault_line = '-W0.5p,black,-'
plot_fault_top_edge = '-W2p,black'
plot_fault_hyp_open = '-W1p,black'


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


################## hf_sim ####################
# check source code in bin dir for more help #
##############################################

# output files in format of: hf_prefix + '_STAT.COMP'
hf_prefix = '22Feb2011_bev01'
# binary that simulates the HF data
hf_sim_bin = '/hpc/home/rwg43/StochSim/Src/V5.4/hb_high_v5.4.4'
# where to save output
hf_outdir = 'Acc'
# duration of HF sim
hf_t_len = '100' # seconds
# HF simulation step. should be small
hf_dt = '0.005' # seconds
# slip model
hf_slip = '/hpc/home/hnr12/RupModel/2011Feb22_m6pt2/Stoch/m6.20-16.0x9.0_s560.stoch'
# 1D velocity model
hf_v_model = '/hpc/home/hnr12/VelocityModel/Mod-1D/Cant1D_v2-midQ.1d'
# for western US, check applicability to NZ
hf_sdrop = '50' # average stress-drop, bars
hf_kappa = '0.045'
hf_fmax = '10'
hf_qfexp = '0.6' # Q freq. exponent
# ???
hf_vs_moho = '999.9'
# uncertainty (sigma) to consider for fourier amplitude
hf_fa_sig_1 = '0.0'
# rupture velocity uncertainty affects corner frequency, scales GMs as f0^2
hf_rv_sig_1 = '0.1'
# ??? 2 rays to consider: direct1 and moho2
hf_rayset = '2 1 2'
# rup vel factor. ratio of rupture : Vs
hf_rvfac = '0.8'
# shallow fault portion multiplier
hf_shal_rvfac = '0.70'
# deep fault portion multiplier
hf_deep_rvfac = '0.70'
# C0 coefficient
hf_czero = '2.1'
# C alpha
hf_calpha = '-99'
# should be '1' (apply BJ97 crustal amplification factors)
# 0 for no site amplification
hf_site_amp = '1'
# seismic moment. '-1': from input rupture model
hf_mom = '-1'
# rupture velocity. '-1': from input rupture model
hf_rupv = '-1.0'
# ??? seed? '0': not random, '1': random
hf_seed = '5481190'




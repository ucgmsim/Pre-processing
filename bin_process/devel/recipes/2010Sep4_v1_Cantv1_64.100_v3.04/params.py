import os.path
from platform import node
from params_base import *

if __name__ == '__main__':
    # in case someone has passed this file as a parameter to python
    print('\nDo not run ' + __file__ + \
            '\nOnly edit variables in ' + __file__ + ' and then run the set param script.\n')
    exit()

# to enable templates, add specific parameters to params_override_<name>.py
# works with binaries/scripts using e3d.par (processed by set_runparams.py)
params_override = '2010Sep4'

# parameters are written to files, need to be strings

######## Global CONSTANTS ########

# things that everyone doesn't have, eg. binaries are within here


# works on Windows and POSIX paths
user_scratch = os.path.join(user_root, 'scratch')


# keep as int, processed before writing to file
# only product is stored
n_proc_x = 8
n_proc_y = 8
n_proc_z = 8

### folowing values are dependent on the velocity model used
# low pass frequency, filter for output seimograms
flo = '1.0' # hertz
# spatial grid spacing
hh = '0.100' # km
# x, y, z grid size (multiples of grid spacing)
nx = '1400'
ny = '1200'
nz = '460'
# model reference location
MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'
###

# cap number of timesteps in simulation, not all timesteps have outputs
# max simulation timeperiod = nt * dt eg: 10,000 * 0.005 = 50 seconds
nt = '20000'
# dt should be 0.005 (or smaller), small increments required in simulation
dt = '0.005'
# how often to save outputs (measured in simulation timesteps)
DUMP_ITINC = '4000' # nt

# output timestep in multiples of simulation timestep
# eg: simulation dt 0.005 sec * dt_ts 20 = 0.1 second increments
dt_ts = '20'
# x, y, z decimation along axis
# store output at lower x, y, z resolution by factor provided
# eg: nx 1400 / dx_ts 5 = 280 points along 1400 nx * 0.1 hh = 140 km
dx_ts = '5'
dy_ts = '5'
dz_ts = '1'

ts_xy='1'
iz_ts='1'

ts_xz='0'
iy_ts='2'

ts_yz='0'
ix_ts='99'


# which time slices to iterate over
ts_start = '0'     # first one is 0
ts_inc = '1'       # increment, larger than 1 to skip
ts_total = '40' # '400'   # number of slices to generate. sim time = ts_total * dt * dt_ts

# swap_bytes 0/1 no/yes - should be 1 if
#   ts_file created on supercomp and this file is run on laptop; zero if run within supercomputer)
# the three lines below not used if the TSFiles are created using 'gen_ts.py' on supercomputer
#   and copied to local computer beforehand.
swap_bytes = '0'
lonlat_out = '1'
scale = '1.0'

# EMOD3D parameters template and complete parameter files
parfile = 'e3d.par' # update load leveler file

#FD_VMODFILE = 'Cant1D_v1.fd_modfile'     #This line was for a 1D Vp,Vs,rho model
# set names of P,S,D files
PMOD = 'vp3dfile.p'
SMOD = 'vs3dfile.s'
DMOD = 'rho3dfile.d'


ENABLE_RESTART = '0'
READ_RESTART = '0'
RESTART_ITINC = '20000'


n_proc = str(n_proc_x * n_proc_y * n_proc_z)

# directories - main. change global_root with user_root as required

stat_dir = os.path.join(global_root, 'StationInfo')
vel_mod_params_dir = os.path.join(global_root, 'VelocityModel/ModelParams')


wcc_prog_dir = os.path.join(global_root, 'EMOD3D/WccFormat/bin/powerpc-AIX-nesi2-xlc')

seis_tmp_dir = os.path.join(user_scratch, run_name, 'SeismoBin')

# directories - derived
#sim_dir = os.path.join(run_dir, run_name)
#hf_sim_dir = os.path.join(run_dir, hf_run_name)
#bb_sim_dir = os.path.join(run_dir, bb_run_name)

restart_dir = os.path.join(lf_sim_dir, 'Restart')
bin_output = os.path.join(lf_sim_dir, 'OutBin')
log_dir = os.path.join(lf_sim_dir, 'Rlog')
slipout_dir = os.path.join(lf_sim_dir,'SlipOut')

# files
stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_nz01-h0.100.statcords')


############# winbin-aio ##############
# -------- START OF USER DEFINED INPUT -------------


# Define location of FD output files to use
SEISDIR = bin_output

# Define folder to write velocity seismograms to
vel_dir = os.path.join(lf_sim_dir,'Vel')

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'

# Define the directory with the station information, and components
FD_STATLIST = os.path.join(stat_dir,'fd_nz01-h0.100.ll')

############### gen_ts ###################

ts_file = os.path.join(bin_output, run_name+ '_xyts.e3d') #the file created by merge_ts

# ABSMAX = 1 all components will be plotted as hypotenuse
#        = 0 plot polar plot of component 0
ABSMAX = '1'

LONLAT_OUT = '1' # LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice

GRIDFILE = os.path.join(vel_mod_params_dir, 'gridout_nz01-h') + hh # GRIDFILE is the file containing the local (x,y,z) coordinates for this 3D run


t_slice_dir = os.path.join(lf_sim_dir, 'TSlice')
ts_out_dir = os.path.join(t_slice_dir, 'TSFiles')
ts_out_prefix = os.path.join(ts_out_dir, run_name)

################## plot_ts ####################

# area to be shown/plotted
# CANTERBURY, SOUTHISLAND
plot_region = 'CANTERBURY'

# main title and subtitles that appears at the top of the picture
plot_main_title = 'Mw7.1 4 Sept 2010 Earthquake'
plot_sub_title = 'Beavan 1Fault, Stoch Slip, v1.64'

# show seismograms for following stations
# longitude, latitude, offset, offset azimuth
# (decimal), (decimal), (km),  (degrees true)
plot_seismo = [ \
[-43.5928, 172.3811, 0, 0], \
[-43.8087, 172.2524, 0, 0], \
[-43.4897, 172.1022, 0, 0]] # ROLC, SBRC, DFHS
# time azimuth (degrees true), x length (km), max y length (km)
plot_seismo_params = [90, 35, 12]
# "SCEC" to make seismograms draw out of origin, "NORM" to make them extend
plot_seismo_style = "SCEC"

# details of the spatial and termporal discretization and spacing
plot_orig_dt = '0.1'  # time step of the time slice output (this is DT*DT_TS from the run files)

# output dirs and resolution (dpi)
plot_ps_dir = os.path.join(t_slice_dir, 'PlotFiles')
plot_png_dir = os.path.join(t_slice_dir, 'Png')
plot_res = '140'   # 720 for PDF quality

# topography
plot_topo_dir = os.path.join(global_root, 'PlottingData/TopoData')
plot_topo_file = os.path.join(plot_topo_dir, 'srtm_71_21.grd')
plot_topo_illu = '-I' + os.path.join(plot_topo_dir, 'srtm_71_21_i5.grd')
plot_topo_file_2 = os.path.join(plot_topo_dir, 'etopo2.grd')
# velocity for plotting (colour palette)
# when absmax=0, max/min is taken as pos/neg
# ie. VALUES HERE SHOULD BE POSITIVE
# min is cutoff (below this is transparent)
plot_topo_a_min = '1.0'
# inc is colour change step
plot_topo_a_inc = '10'
# anything at or above max shows up as max
plot_topo_a_max = '50'

# 'finitefault' or 'pointsource'
plot_type = 'finitefault'
# fault plane resources
fault_plane_dir = os.path.join(global_root, 'PlottingData', 'sourcesAndStrongMotionStations')
fault_file = os.path.join(fault_plane_dir, 'bev01_DarfieldFaultPlane.xy')
plot_fault_add_plane = os.path.join(fault_plane_dir, 'addStandardFaultPlane.sh')
plot_fault_line = '-W0.5p,black,-'
plot_fault_top_edge = '-W2p,black'
plot_fault_hyp_open = '-W1p,black'
# pointsource resources
plot_beachball = '171.9773 -43.252 8 154 83 16 62 74 173 7.98 22 171.9773 -43.252 4.6'

# set PALETTE = '-Crelief.cpt'
plot_palette = '-Cgray'

# location of offset plotting (for when 3 component plotting used - not currently utilized)
plot_x_org = '1.15'
plot_y_org = '2.5'
plot_x_inch = '5.0'
plot_x_shift = '%0.6f' % (float(plot_x_inch) + 0.2)

# specify the increments of X/Y (cartesian coords) for masks etc.
plot_dx = '0.002'
plot_dy = '0.002'


################## hf_sim ####################
# check source code in bin dir for more help #
##############################################

# binary that simulates the HF data
hf_sim_bin = os.path.join(global_root,'EMOD3D/StochSim/Src/V5.4/hb_high_v5.4.5')
hf_prefix = 'hf'
# duration of HF sim
hf_t_len = str(int(int(nt) * float(dt))) # seconds
# HF simulation step. should be small
hf_dt = dt # seconds
# slip model
# 1D velocity model
# for western US, check applicability to NZ
hf_sdrop = '50' # average stress-drop, bars
hf_kappa = '0.045'
hf_fmax = '10'
hf_qfexp = '0.6' # Q freq. exponent
# depth to moho
hf_vs_moho = '999.9'
# uncertainty (sigma) to consider for fourier amplitude
hf_fa_sig_1 = '0.0'
# rupture velocity uncertainty affects corner frequency, scales GMs as f0^2
hf_rv_sig_1 = '0.1'
# 2 rays to consider: direct1 and moho2
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
# seed '0': not random, '1': random
hf_seed = '5481190'


############ acc2vel and match_seismo #############

# binary that integrates / differentiates (Acc <> Vel)
int_bin = os.path.join(wcc_prog_dir,'integ_diff')

# binary that applies site amplification
siteamp_bin = os.path.join(wcc_prog_dir,'wcc_siteamp')

# binary that applies ??? filter
tfilter_bin = os.path.join(wcc_prog_dir,'wcc_tfilter')
# binary that adds LF and HF using matched filters
match_bin = os.path.join(wcc_prog_dir,'wcc_add')

# binary which returns PGA
getpeak_bin = os.path.join(wcc_prog_dir,'wcc_getpeak')
# components which will be converted (acc2vel)
int_comps = ['000', '090', 'ver']
match_log = 'seismo.log'
# station file with vs30 reference values
stat_vs_ref = os.path.join(stat_dir,'cantstations_cant1D_v1.vs30ref')
# station file with vs30 estimates
stat_vs_est = os.path.join(stat_dir,'cantstations.vs30')
# filter properties for HF (short period) ground motion
match_hf_fhi = '1.0'        # high pass frequency, Hz
match_hf_flo = '1.0e+15'    # low pass frequency, Hz
match_hf_ord = '4'          # filter order
match_hf_tstart = '0.0'     # for alignment
# filter properties for LF (long period) ground motion
match_lf_fhi = '0.0'        # high pass frequency, Hz
match_lf_flo = '1.0'        # low pass frequency, Hz
# '4': standard?, '0': already low pass filtered
match_lf_ord = '0'          # filter order
match_lf_tstart = '0.0'     # for alignment
# list of matching HF and LF components
match_hf_comps = ['000', '090', 'ver']
match_lf_comps = ['000', '090', 'ver']
# vs30 amplification model [cb2008 | bssa2014 | cb2014]
site_amp_model = 'cb2014'
# relating to site amplification
site_vref_max = '1100'  # reference vs30 for cb08/14 amp model
site_fmin = '0.1'       # 0.2 Hz = 5 sec. point where tapering to unity begins
site_fmidbot = '0.2'    # freq. for which cap is applied f = 1 Hz => T = 1 sec in GP10
site_flowcap = '0.0'
# set to vs30 used in 1D model for high frequency runs (VREF for HF)
GEN_ROCK_VS = 500


########### gen_statgrid ##########

# any values to omit near the boundaries
X_BND_PAD = 0
Y_BND_PAD = 0
# input to process
MODEL_COORDS = os.path.join(vel_mod_params_dir, 'model_coords_nz01-h0.100')
# output statgrid file
STATGRID_GEN = 'statgrid-0.5x0.5-nz01_h0.100.ll'





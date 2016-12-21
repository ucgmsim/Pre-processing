import os.path
from platform import node
from params_base import *

import params_vel as vel 
import params_src as src
import params_stat as stat

# number of processes
n_proc = '512'  
# low pass frequency, filter for output seimograms
flo = '1.0' # hertz  

######## Global CONSTANTS ########
#FROM RUPTURE MODEL
# model reference location
MODEL_LAT = vel.ORIGIN_LAT
MODEL_LON = vel.ORIGIN_LON

# dt should be 0.005 (or smaller), small increments required in simulation
dt = src.DT

#FROM VELOCITY MODEL
MODEL_ROT = vel.ORIGIN_ROT
# spatial grid spacing
hh = vel.HH
# x, y, z grid size (multiples of grid spacing)
nx = vel.NX
ny = vel.NY
nz = vel.NZ

sufx = vel.SUFX #need to get from VELOCITY MODEL

# cap number of timesteps in simulation, not all timesteps have outputs
# max simulation timeperiod = nt * dt eg: 10,000 * 0.005 = 50 seconds
sim_duration = vel.SIM_DURATION



#############################################################################
# works on Windows and POSIX paths
user_scratch = os.path.join(user_root, 'scratch')

nt = str(int(float(sim_duration)/float(dt)))
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
ts_inc = '1'       # increment, larger than 1 to skip
ts_start = '0'
ts_total = str(int(float(sim_duration)/(float(dt)*float(dt_ts))))

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


# directories - main. change global_root with user_root as required

wcc_prog_dir = os.path.join(global_root, 'EMOD3D/WccFormat/bin/powerpc-AIX-nesi2-xlc')
seis_tmp_dir = os.path.join(user_scratch, run_name, 'SeismoBin')

# directories - derived
#sim_dir = os.path.join(run_dir, run_name)
#hf_sim_dir = os.path.join(run_dir, hf_run_name)
#bb_sim_dir = os.path.join(run_dir, bb_run_name)

#restart_dir = os.path.join(lf_sim_dir, 'Restart') #see params_uncertain.py
#bin_output = os.path.join(lf_sim_dir, 'OutBin') #see params_uncertain.py
#log_dir = os.path.join(lf_sim_dir, 'Rlog') #see params_uncertain.py
#slipout_dir = os.path.join(lf_sim_dir,'SlipOut') #see params_uncertain.py


############# winbin-aio ##############
# -------- START OF USER DEFINED INPUT -------------

# Define location of FD output files to use
#SEISDIR = bin_output #see params_uncertain.py

# Define folder to write velocity seismograms to
# vel_dir = os.path.join(lf_sim_dir,'Vel') #see params_uncertain.py

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'


############### gen_ts ###################

#ts_file = os.path.join(bin_output, run_name+ '_xyts.e3d') #the file created by merge_ts #see params_uncertain.py

# TODO: not used anywhere???
# ABSMAX =1 vector magnitude of all 3 components is output into <fileroot>.0
#        =0 3 components are output respectively into <fileroot>.0, <fileroot>.1, <fileroot>.2
ABSMAX = '1'

LONLAT_OUT = '1' # LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice


#t_slice_dir = os.path.join(lf_sim_dir, 'TSlice') #see params_uncertain.py
#ts_out_dir = os.path.join(t_slice_dir, 'TSFiles') #see params_uncertain.py
#ts_out_prefix = os.path.join(ts_out_dir, run_name) #see params_uncertain.py



################## hf_sim ####################
# check source code in bin dir for more help #
##############################################

# binary that simulates the HF data
hf_sim_bin = os.path.join(global_root,'EMOD3D/StochSim/Src/V5.4/hb_high_v5.4.5_np2mm+')
hf_prefix = 'hf'
# duration of HF sim
hf_t_len = '100' # seconds
# HF simulation step. should be small
hf_dt = '0.005' # seconds
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
site_fmin = '0.2'       # 0.2 Hz = 5 sec. point where tapering to unity begins
site_fmidbot = '0.5'    # freq. for which cap is applied f = 1 Hz => T = 1 sec in GP10
site_flowcap = '0.0'
site_fmid = '1.0'
site_fhigh = '3.333'
site_fhightop = '10.0'
site_fmax = '15.0'

# set to vs30 used in 1D model for high frequency runs (VREF for HF)
GEN_ROCK_VS = 500


########### gen_statgrid ##########

# any values to omit near boundaries
# used when generating the statgrid
X_BND_PAD = 0
Y_BND_PAD = 0
# stations to omit near boundaries
# used when processing stations > grid coords
X_BND = 0 # original was 30
Y_BND = 0 # original was 30
# wcc options
CENTER_ORIGIN = '1'
# distance calculation
GEOPROJ = '1'
# binaries
XYZ2LL_BIN = '/nesi/projects/nesi00213/tools/xyz2lonlat'
LL2GRID_BIN = '/nesi/projects/nesi00213/tools/latlon2statgrid'
# station file list for generating station grid
STAT_FILES = [stat_file]

STATGRID_GEN = 'statgrid-%sx%s-%s' %(float(hh)*float(dx_ts),float(hh)*float(dy_ts),sufx) # output statgrid file


########### stat_max.py #########

# where maximum values are stored
MAX_STATFILE = 'stations_max.bin'


if __name__ == '__main__':
    # in case someone has passed this file as a parameter to python
    print('\nDo not run ' + __file__ + \
            '\nOnly edit variables in ' + __file__ + ' and then run the set param script.\n')
    exit()


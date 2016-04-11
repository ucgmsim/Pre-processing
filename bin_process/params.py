import os.path
# parameters are written to files, need to be strings

######## Global CONSTANTS ########

RUN_NAME = '2011Feb22_m6pt2bev01_Cantv1.64'
VERSION = '3.0.4'
EXTENDED_RUN_NAME = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_Test'


NPROC_X = 16
NPROC_Y = 16
NPROC_Z = 8


FLO = '1.0'
HH = '0.100'
NX = '1400'
NY = '1200'
NZ = '460'
NT = '20000'
DT = '0.005'

MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'

DUMP_ITINC = '4000'

DT_TS = '20'
DX_TS = '5'
DY_TS = '5'
DZ_TS = '1'


DEFAULT_PARFILE = 'e3d_default.par'
PARFILE = 'e3d.par'

#FD_VMODFILE = 'Cant1D_v1.fd_modfile'     #This line was for a 1D Vp,Vs,rho model
# set names of P,S,D files
PMOD = 'vp3dfile.p'
SMOD = 'vs3dfile.s'
DMOD = 'rho3dfile.d'

VMOD_VER = 'v1.64'

ENABLE_RESTART = '0'
READ_RESTART = '0'
RESTART_ITINC = '20000'


######## Variables derived from Global CONSTANTS ##########

# where is this used? needs to be a string?
NPROC = NPROC_X * NPROC_Y * NPROC_Z


ROOT = '/hpc/scratch/nesi00213'
RUN_DIR_ROOT = ROOT + '/RunFolder'
SRF_DIR_ROOT = ROOT + '/RupModel'
STAT_DIR = ROOT + '/StationInfo'
MOD_DIR_ROOT = ROOT


# XXX: was 3.04 while version elsewhere was 3.0.4, changing this changed a folder
#        while changing version changes parameters

SIMDIR_ROOT = RUN_DIR_ROOT + '/'+EXTENDED_RUN_NAME
MAIN_OUTPDIR = SIMDIR_ROOT + '/OutBin'

VMODDIR_ROOT = ROOT+'/CanterburyVelocityModel'

SRF_FILE = SRF_DIR_ROOT + '/Srf/m6.20-16.0x9.0_s560.srf' #rmc
STATCORDS = STAT_DIR + '/fd_nz01-h0.100.statcords'

#WCC_PROGDIR = ROOT + '/EMOD3D/WccFormat/bin/powerpc-AIX-p1n14-c-xlc'
WCC_PROGDIR = ROOT + '/EMOD3D/WccFormat/bin/ppc64-Linux-p2n14-c-gcc'

MAIN_RESTARTDIR = SIMDIR_ROOT + '/Restart'


VMODDIR = VMODDIR_ROOT +'/'+ VMOD_VER
VELMODPARAMSDIR = ROOT + '/VelocityModel/ModelParams'

TMP_SEISDIR = '/hpc/scratch/' + os.getenv('USER') + SIMDIR_ROOT + '/SeismoBin'




############# winbin-aio ##############
# -------- START OF USER DEFINED INPUT -------------


# Define location of FD output files to use
SEISDIR = MAIN_OUTPDIR 

# Define folder to write velocity seismograms to
VELDIR = os.path.join(SIMDIR_ROOT,'Vel')

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'
SCALE = '1.0'

# Define the directory with the station information, and components
FD_STATLIST = STAT_DIR + '/fd_nz01-h0.100.ll'


############### gen_ts ###################

SIMDIR = MAIN_OUTPDIR #the directory name where the 3D FD run output exists
TSFILE = os.path.join(SIMDIR, EXTENDED_RUN_NAME + '_xyts.e3d') #the file created by merge_ts

HH = '0.100' # HH is the grid spacing used for the 3D FD run
DXTS = '5' # DXTS is the decimation along the X-axis for the time slice output (same as dxts in e3d.par)
DYTS = '5' # DYTS is the decimation along the Y-axis for the time slice output
DZTS = '1' # DZTS is the decimation along the Z-axis for the time slice output (=1 for xy time slice)


# ABSMAX =1 vector magnitude of all 3 components is output into <fileroot>.0
#        =0 3 components are output respectively into <fileroot>.0, <fileroot>.1, <fileroot>.2
ABSMAX = '1'

TS_INC = 1 # TS_INC is the increment for output time slices (=2 means every other slice is output)
TS_START = 0 # TS_START is the starting time slice (0 is first one)
TS_TOTAL = 400 # the total number of slices to generate. the total sim time = TS_TOTAL * ORIG_DT


LONLAT_OUT = '1' # LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice
GRIDFILE = VELMODPARAMSDIR + '/gridout_nz01-h' + HH # GRIDFILE is the file containing the local (x,y,z) coordinates for this 3D run

TSLICE_DIR = SIMDIR_ROOT+'/TSlice'
TS_OUTFILE_DIR = os.path.join(TSLICE_DIR,'TSFiles')
TS_OUTFILE_PREFIX = os.path.join(TS_OUTFILE_DIR, RUN_NAME)

SWAP_BYTES = '0'
SCALE = '1.0'


################## plot_ts ####################


# main title and subtitles that appears at the top of the picture
MAIN_TITLE = 'Mw6.2 22 Feb 2011 Earthquake'
SUB_TITLE = 'Beavan 1Fault, Stoch Slip, Chch 1D VM'


## high-level directories
#MAINDIR = os.path.expanduser('~bab70')
## not used?
#SIMDIR = '../'

# input can be provided either via the "._xyts.e3d" file (option=1)
# or via all of the individual TSFiles, which are obtained from gen_ts.csh 
# within BlueFern and then scp'd to the local machine (option=2)
option = 2

# information required for option = 1
# ----------------------------------


#
DXTS = '5'  #reduce X and Y resolution by using these factors
DYTS = '5'

# information required for option 2.
# ----------------------------------

# details of the spatial and termporal discretization and spacing
HH = '0.100'  # spatial grid spacing in km
ORIG_DT = '0.1'  # time step of the time slice output (this is DT*DT_TS from the run files)
MODELPARAMS = VELMODPARAMSDIR + '/model_params_nz01-h' + HH  # used to get location of model corners


# specify whether or not to byte-swap (=0 no; =1 yes - which should be used if
#   TSFILE created on supercomp and this file is run on laptop; zero if run within supercomputer)
# the three lines below not used if the TSFiles are created using 'gen_ts.py' on supercomputer 
#   and copied to local computer beforehand.
SWAP_BYTES = '0'
LONLAT_OUT = '1'
SCALE = '1.0'

#components to consider (0= when only looking at vector magnitude [i.e. ABSMAX=1 further down].
COMPS = [ 0 ]


#specify the directory for PNG files and the resolution (dpi) to create
PSDIR = TSLICE_DIR  + '/PlotFiles'
PNGDIR = TSLICE_DIR  + '/Png'
RES = '140'   #720 for PDF quality
# make directories for storing the .ps and .png files
if not os.path.exists(PSDIR):
    os.makedirs(PSDIR)
if not os.path.exists(PNGDIR):
    os.makedirs(PNGDIR)

TOPODIR = ROOT + '/PlottingData/TopoData'
TOPO_FILE = TOPODIR + '/srtm_71_21.grd'
ILLU = '-I' + TOPODIR + '/srtm_71_21_i5.grd'
TOPO_FILE2 = TOPODIR + '/etopo2.grd'

# set PALETTE = '-Crelief.cpt'
PALETTE = '-Cgray'

# specific locations to display on figure (not currently used)
SITES = [ 'Rolleston', 'Darfield', 'Lyttelton', 'Akaroa', 'Kaiapoi', 'Rakaia', 'Oxford' ]
# alignment; 2letter, L,C,R (left, center, right); T,M,B (top, middle, bottom)
SPOS = [ 'RB', 'CB', 'LM', 'RB', 'LB', 'RT', 'LB' ]
SLON = [ '172.3791667', '172.1116667', '172.7194444', '172.9683333',
        '172.6569444', '172.0230556', '172.1938889' ]
SLAT = [ '-43.59083333', '-43.48972222', '-43.60305556', '-43.80361111',
        '-43.38277778', '-43.75611111', '-43.29555556' ]
# specifying plotting preferences for site locations
SSYM = 'c0.10' # symbol
SFIL = '220/220/220' # fill color
SLIN = '1,000/000/000' # line?

# location of offset plotting (for when 3 component plotting used - not currently utilized)
XORG = '1.15'
YORG = '2.5'
XINCH = 5.0
XSHFT = '%0.6f' % (XINCH + 0.2)

# specify the maximum plotting window for the basemap
PLOT_XMIN = '171.75'
PLOT_XMAX = '173.00'
PLOT_YMIN = '-44.00'
PLOT_YMAX = '-43.20'
# region in GMT format
PLOT_REGION = '/'.join([PLOT_XMIN, PLOT_XMAX, PLOT_YMIN, PLOT_YMAX])

# specify the plotting region for the time slice

TS_XMIN = '171.75'
TS_XMAX = '173.00'
TS_YMIN = '-44.00'
TS_YMAX = '-43.20'
# region in GMT format
TS_REGION = '/'.join([TS_XMIN, TS_XMAX, TS_YMIN, TS_YMAX])
# specify the increments of X/Y (cartesian coords) for masks etc.
DX = '0.002'
DY = '0.002'


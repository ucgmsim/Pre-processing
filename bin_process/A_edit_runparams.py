# not executable #!/usr/bin/env python2
# edit all parameters from this file
# warning: this filename is referenced by the set parameters script

import os

if __name__ == '__main__':
    # in case someone has passed this file as a parameter to python
    print('\nDo not run ' + __file__ + \
            '\nOnly edit variables in ' + __file__ + ' and then run the set param script.\n')
    exit()

# ================================================================ #
#                     EDIT FOLLOWING VALUES                        #
# ================================================================ #
#       parameters have to be strings (enclosed in '' or "")       #

VERSION = '3.0.4'

# currently not used?
NPROC_X = 16
NPROC_Y = 16
NPROC_Z = 8
# where is this used? should update ll file?
NPROC = NPROC_X * NPROC_Y * NPROC_Z

FLO = '1.0'
HH = '0.100'
# x, y, z distance (km)
NX = '1400'
NY = '1200'
NZ = '460'
# number of timesteps n, time step Dt
NT = '20000'
DT = '0.005'

RUN_NAME = '2011Feb22_m6pt2bev01_Cantv1.64'
RUN_DIR_ROOT = os.path.expanduser('~rmc84') + '/RunFolder'
SRF_DIR_ROOT = os.path.expanduser('~rmc84') + '/RupModel'
MOD_DIR_ROOT = os.path.expanduser('~rmc84')

# XXX: was 3.04 while version elsewhere was 3.0.4, changing this changed a folder
#        while changing version changes parameters
SIMDIR_ROOT = RUN_DIR_ROOT + '/LPSim-2011Feb22b560_v1_Cantv1_64-h0.100_v3.04'
MAIN_OUTPDIR = SIMDIR_ROOT + '/OutBin'
VMODDIR_ROOT = os.path.expanduser('~emt41') + '/CanterburyVelocityModel'
SRF_FILE = SRF_DIR_ROOT + '/Srf/m6.20-16.0x9.0_s560.srf' #rmc
STATCORDS = MOD_DIR_ROOT + '/StationInfo/fd_nz01-h0.100.statcords'

# using default as template, additional parameters are then appended
DEFAULT_PARFILE = 'e3d_default.par'
PARFILE = 'e3d.par'

MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'

DUMP_ITINC = '4000'

DT_TS = '20'
DX_TS = '5'
DY_TS = '5'
DZ_TS = '1'

ENABLE_RESTART = '0'
READ_RESTART = '0'
MAIN_RESTARTDIR = SIMDIR_ROOT + '/Restart'
RESTART_ITINC = '20000'

#FD_VMODFILE = 'Cant1D_v1.fd_modfile'     #This line was for a 1D Vp,Vs,rho model
# set names of P,S,D files
PMOD = 'vp3dfile.p'
SMOD = 'vs3dfile.s'
DMOD = 'rho3dfile.d'


SEISDIR = 'SeismoBin'
VMODDIR = VMODDIR_ROOT + '/v1.64'
TMP_SEISDIR = '/hpc/scratch/' + os.getenv('USER') + '/' + SIMDIR_ROOT + '/' + SEISDIR


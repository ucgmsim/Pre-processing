#!/usr/bin/env python

import fileinput
import os
import sys
from __future__ import print_function
from shutil import copyfile

def change_var(filename, variable, value):
    for line in fileinput.input(filename, inplace = True):
        if line.startswith(variable + ' = '):
            line = variable + ' = ' + value + '\n'
        # don't include extra '\n'
        sys.stdout.write(line)

# parameters are written to files, need to be strings
VERSION = '3.0.4'

NPROC_X = 16
NPROC_Y = 16
NPROC_Z = 8
# where is this used? needs to be a string?
NPROC = NPROC_X * NPROC_Y * NPROC_Z

FLO = '1.0'
HH = '0.100'
NX = '1400'
NY = '1200'
NZ = '460'
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
VMODDIR_ROOT = '/hpc/home/emt41/CanterburyVelocityModel'
SRF_FILE = SRF_DIR_ROOT + '/Srf/m6.20-16.0x9.0_s560.srf' #rmc
STATCORDS = MOD_DIR_ROOT + '/StationInfo/fd_nz01-h0.100.statcords'

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
TMP_SEISDIR = '/hpc/scratch/' + os.genenv('USER') + '/' + SIMDIR_ROOT + '/' + SEISDIR

try:
    copyfile(DEFAULT_PARFILE, PARFILE)
except IOError:
    print('Cannot copy DEFAULT_PARFILE to PARFILE! ' + __file__, file = sys.stderr)
    raise

try:
    par_handle = open(PARFILE, 'a')
except IOError:
    print('PARFILE cannot be opened to append data' + __file__, file = sys.stderr)
    raise

par_handle.write('version=' + VERSION + '-mpi')
par_handle.write('name=' + RUN_NAME)
par_handle.write('nx=' + NX)
par_handle.write('ny=' + NY)
par_handle.write('nz=' + NZ)
par_handle.write('h=' + HH)
par_handle.write('nt=' + NT)
par_handle.write('dt=' + DT)

par_handle.write('bfilt=4')
par_handle.write('flo=' + FLO)
par_handle.write('fhi=0.0')
par_handle.write('bforce=0')
par_handle.write('pointmt=0')
par_handle.write('dblcpl=0')
par_handle.write('ffault=2')
par_handle.write('faultfile=' + SRF_FILE)

par_handle.write('model_style=1')
# only for the 1D velocity model
#par_handle.write('model=' + FD_VMODFILE)
par_handle.write('vmoddir=' + VMODDIR)
par_handle.write('pmodfile=' + PMOD)
par_handle.write('smodfile=' + SMOD)
par_handle.write('dmodfile=' + DMOD)
par_handle.write('qpfrac=100')
par_handle.write('qsfrac=50')
par_handle.write('qpqs_factor=2.0')
par_handle.write('fmax=25.0')
par_handle.write('fmin=0.01')
par_handle.write('vmodel_swapb=1')

par_handle.write('modellon=' + MODEL_LON)
par_handle.write('modellat=' + MODEL_LAT)
par_handle.write('modelrot=' + MODEL_ROT)

par_handle.write('enable_output_dump=1')
par_handle.write('dump_itinc=' + DUMP_ITINC)
par_handle.write('main_dump_dir=' + MAIN_OUTPDIR)
par_handle.write('nseis=1')
par_handle.write('seiscords=' + STATCORDS)
par_handle.write('seisdir=' + TMP_SEISDIR)

par_handle.write('ts_xy=1')
par_handle.write('iz_ts=1')
par_handle.write('dtts=' + DT_TS)
par_handle.write('dxts=' + DX_TS)
par_handle.write('dyts=' + DY_TS)
par_handle.write('dzts=' + DZ_TS)

par_handle.write('enable_restart=' + ENABLE_RESTART)
par_handle.write('restartdir=' + MAIN_RESTARTDIR)
par_handle.write('restart_itinc=' + RESTART_ITINC)
par_handle.write('read_restart=' + READ_RESTART)
par_handle.write('restartname=' + RUN_NAME)

par_handle.close()

############# updating variables in other files #############
change_var('winbin-aio.csh', 'RUN', RUN_NAME)
change_var('merge_tsP3.py', 'FILEROOT', RUN_NAME)
change_var('gen_ts_default.py', 'FILEROOT', RUN_NAME)
change_var('plot_ts_bluefern.py', 'NAME', RUN_NAME)


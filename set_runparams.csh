#!/bin/csh

set VERSION = 3.0.4

set NPROC_X = 36
set NPROC_Y = 14
set NPROC_Z = 4
set NPROC = `echo $NPROC_X $NPROC_Y $NPROC_Z | gawk '{printf "%.0f\n",$1*$2*$3;}'`

set FLO = 0.25
set HH = 0.400
set NX = 2125
set NY = 750
set NZ = 250
set NT = 60000
set DT = 0.005

set RUN_NAME = AlpineM7.9SthHypo_v1_SthIslandv1.01-h0.400_v3.04
set RUN_DIR_ROOT = /hpc/home/bab70/RunFolder
set SRF_DIR_ROOT = /hpc/home/bab70/RupModel
set MOD_DIR_ROOT = /hpc/home/bab70

set SIMDIR_ROOT = ${RUN_DIR_ROOT}/LFSim-AlpineM7.9SthHypo_v1_SthIslandv1.01-h0.400_v3.04
set MAIN_OUTPDIR = ${SIMDIR_ROOT}/OutBin
set VMODDIR_ROOT = ${MOD_DIR_ROOT}/VelocityModel
set SRF_FILE = ${SRF_DIR_ROOT}/AlpineFaultTest/Srf/m7.90-411.0x17.3_s1129571.srf
set STATCORDS = ${MOD_DIR_ROOT}/StationInfo/fd_sinz01-h0.400.statcords

set DEFAULT_PARFILE = e3d_default.par
set PARFILE = e3d.par

set MODEL_LAT = -43.75
set MODEL_LON = 170.75
set MODEL_ROT = -50.0

set DUMP_ITINC = 4000

set DT_TS = 20
set DX_TS = 5
set DY_TS = 5
set DZ_TS = 1

set ENABLE_RESTART = 0
set READ_RESTART = 0
set MAIN_RESTARTDIR = ${SIMDIR_ROOT}/Restart
set RESTART_ITINC = 20000

set FD_VMODFILE = Cant1D_v1.fd_modfile

set SEISDIR = SeismoBin
set VMODDIR = ${VMODDIR_ROOT}/Mod-1D
set TMP_SEISDIR = /hpc/scratch/${USER}/${SIMDIR_ROOT}/${SEISDIR}

\cp $DEFAULT_PARFILE $PARFILE

echo "version=${VERSION}-mpi"           >> $PARFILE
echo "name=${RUN_NAME}"			>> $PARFILE
echo "nx=$NX"                           >> $PARFILE
echo "ny=$NY"                           >> $PARFILE
echo "nz=$NZ"                           >> $PARFILE
echo "h=$HH"                            >> $PARFILE
echo "nt=$NT"                           >> $PARFILE
echo "dt=$DT"                           >> $PARFILE

echo "nproc_x=$NPROC_X"                           >> $PARFILE
echo "nproc_y=$NPROC_Y"                           >> $PARFILE
echo "nproc_z=$NPROC_Z"                           >> $PARFILE

echo "bfilt=4"				>> $PARFILE
echo "flo=$FLO"				>> $PARFILE
echo "fhi=0.0"				>> $PARFILE
echo "bforce=0"				>> $PARFILE
echo "pointmt=0"			>> $PARFILE
echo "dblcpl=0"                         >> $PARFILE
echo "ffault=2"                         >> $PARFILE
echo "faultfile=${SRF_FILE}"            >> $PARFILE

echo "model_style=0"                    >> $PARFILE
echo "model=$FD_VMODFILE"               >> $PARFILE
echo "vmoddir=${VMODDIR}"               >> $PARFILE
echo "qpfrac=-1"                        >> $PARFILE
echo "qsfrac=-1"                        >> $PARFILE
echo "qpqs_factor=2.0"                  >> $PARFILE
echo "fmax=25.0"				>> $PARFILE
echo "fmin=0.01"			>> $PARFILE

echo "modellon=$MODEL_LON"              >> $PARFILE
echo "modellat=$MODEL_LAT"              >> $PARFILE
echo "modelrot=$MODEL_ROT"              >> $PARFILE

echo "enable_output_dump=1"		>> $PARFILE
echo "dump_itinc=$DUMP_ITINC"		>> $PARFILE
echo "main_dump_dir=${MAIN_OUTPDIR}"    >> $PARFILE
echo "nseis=1"				>> $PARFILE
echo "seiscords=${STATCORDS}"           >> $PARFILE
echo "seisdir=${TMP_SEISDIR}"           >> $PARFILE

echo "ts_xy=1"				>> $PARFILE
echo "iz_ts=1"				>> $PARFILE
echo "dtts=$DT_TS"			>> $PARFILE
echo "dxts=$DX_TS"			>> $PARFILE
echo "dyts=$DY_TS"			>> $PARFILE
echo "dzts=$DZ_TS"			>> $PARFILE

echo "enable_restart=${ENABLE_RESTART}" >> $PARFILE
echo "restartdir=${MAIN_RESTARTDIR}"    >> $PARFILE
echo "restart_itinc=${RESTART_ITINC}"   >> $PARFILE
echo "read_restart=${READ_RESTART}"     >> $PARFILE
echo "restartname=${RUN_NAME}"          >> $PARFILE

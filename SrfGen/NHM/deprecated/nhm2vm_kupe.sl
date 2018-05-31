#!/bin/bash
# script version: slurm

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=nhm2vm
#SBATCH --account=nesi00213
#SBATCH --partition=NeSI
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output nhm2vm-%j.out
#SBATCH --error nhm2vm-%j.err
###SBATCH --mail-type=all
###SBATCH --mail-user=test@test.com
###SBATCH --mem-per-cpu=16G
###SBATCH -C avx
#SBATCH --hint=nomultithread

## END HEADER
date

source machine_env.sh
export SRFGENPATH=/nesi/nobackup/nesi00213/tmp/Pre-processing/SrfGen
export PYTHONPATH=/nesi/nobackup/nesi00213/tmp/qcore:/$PYTHONPATH:$SRFGENPATH:/nesi/nobackup/nesi00213/tmp/post-processing/im_processing/computations/GMPE_models
export RUNDIR=/nesi/nobackup/nesi00213/tmp/auto_preproc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nesi/nobackup/nesi00213/gmt/gmt-stable/lib64
export HDF5_USE_FILE_LOCKING=FALSE # to fix HDF5 file locking issue https://github.com/ALPSCore/ALPSCore/issues/348
#export PMI_MMAP_SYNC_WAIT_TIME=300

srun python2 $SRFGENPATH/NHM/nhm2vm_mp.py $RUNDIR/nhm_selection_sample.txt
#python2 $SRFGENPATH/NHM/nhm2vm_mp.py $RUNDIR/nhm_selection_geo_si.txt

date


#!/bin/bash
# script version: slurm

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=nhm2srf
#SBATCH --account=nesi00213
#SBATCH --partition=NeSI
#SBATCH --ntasks=40
#SBATCH --time=00:20:00
#SBATCH --output nhm2srf.out
#SBATCH --error nhm2srf.err
###SBATCH --mail-type=all
###SBATCH --mail-user=test@test.com
###SBATCH --mem-per-cpu=16G
###SBATCH -C avx
#SBATCH --hint=nomultithread

## END HEADER
source machine_env.sh
date
export SRFGENPATH=/nesi/nobackup/nesi00213/Pre-processing/SrfGen
export PYTHONPATH=$PYTHONPATH:$SRFGENPATH
export RUNDIR=/nesi/nobackup/nesi00213/tmp/auto_preproc
#need to have a copy of NZ_FLTmodel_2010.txt in $RUNDIR
srun python2 $SRFGENPATH/NHM/nhm2srf_hpc.py $RUNDIR/nhm_selection_geo_si.txt  
date


#!/bin/bash
# script version: slurm

# sbatch NHM2SRF_VM.sl /path/to/NHM_SELECTION_FILE [-p]
# 

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=nhm2srf_vm
#SBATCH --account=nesi00213
#SBATCH --partition=NeSI
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --output nhm2srf_vm-%j.out
#SBATCH --error nhm2srf_vm-%j.err
###SBATCH -C avx
#OpenMP+Hyperthreading works well for VM

## END HEADER
source machine_env.sh

# user specific parameters
export SRF_NUMPROCS=72
export VM_NUMPROCS=4
export OMP_NUM_THREADS=18

export SELECTION_FILE=$1
# important paths
export SRFGENPATH=/nesi/project/nesi00213/Pre-processing/SrfGen

# leave these alone
export PYTHONPATH=$PYTHONPATH:$SRFGENPATH

 #to fix HDF5 file locking issue https://github.com/ALPSCore/ALPSCore/issues/348
export HDF5_USE_FILE_LOCKING=FALSE

echo "SRF generation starts"
date
python $SRFGENPATH/NHM/nhm2srf.py $SELECTION_FILE -n $SRF_NUMPROCS
echo "SRF generation completed"
date
echo "VM generation starts"
srun python $SRFGENPATH/srfinfo2vm.py "autosrf/*/Srf/*.info" --pgv 2 --hh 0.4 --dt 0.02 --space-land 5 --space-srf 15 --min-vs 0.5 -n $VM_NUMPROCS
date
echo "VM generation completed"


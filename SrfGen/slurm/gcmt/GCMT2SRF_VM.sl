#!/bin/bash
# script version: slurm

# sbatch GCMT2SRF_VM.sl /path/to/GCMT_FILE [/path/to/fault_file /path/to/uncertainty.yaml]

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=gcmt2srf_vm
#SBATCH --account=nesi00213
#SBATCH --partition=large
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output gcmt2srf_vm-%j.out
#SBATCH --error gcmt2srf_vm-%j.err
# OpenMP+Hyperthreading works well for VM
###SBATCH --hint=nomultithread

## END HEADER

# user specific parameters
export SRF_NUMPROCS=72
export VM_NUMPROCS=4
export OMP_NUM_THREADS=18

export GCMT_FILE=$1

# important paths
export SRFGENPATH=/nesi/project/nesi00213/Pre-processing/SrfGen

if [ ! -z $3 ]; then uncertainty="-r $3"; fi
if [ ! -z $2 ]; then fault_file="-c $2"; fi

 #leave these alone
export PYTHONPATH=$PYTHONPATH:$SRFGENPATH
# to fix HDF5 file locking issue https://github.com/ALPSCore/ALPSCore/issues/348
export HDF5_USE_FILE_LOCKING=FALSE

echo "SRF generation starts"
date
python $SRFGENPATH/gcmt2srf.py $GCMT_FILE $SRFGENPATH/lp_generic1d-gp01_v1.vmod -n $SRF_NUMPROCS $uncertainity $fault_file
echo "SRF generation completed"
date
echo "VM generation starts"
srun python $SRFGENPATH/srfinfo2vm.py "autosrf/*/Srf/*.info" -n $VM_NUMPROCS --pgv -1.0 --hh 0.4 --dt 0.02 --min-vs 0.5 --vm-version "1.66" --vm-topo "SQUASHED_TAPERED"
date
echo "VM generation completed"

## Only one VM is generated 3631755 as the other one (2122842) is in the ocean 

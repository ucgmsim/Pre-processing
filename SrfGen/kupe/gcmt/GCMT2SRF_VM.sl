#!/bin/bash
# script version: slurm

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=gcmt2srf_vm
#SBATCH --account=nesi00213
#SBATCH --partition=NeSI
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output gcmt2srf_vm-%j.out
#SBATCH --error gcmt2srf_vm-%j.err
###SBATCH --mail-type=all
###SBATCH --mail-user=test@test.com
###SBATCH --mem-per-cpu=16G
###SBATCH -C avx
#SBATCH --hint=nomultithread

## END HEADER
source machine_env.sh

# user specific parameters
export NUMPROCS=40
export GCMT_FILE=`pwd`/GeoNet_CMT_solutions.csv

# important paths
export SRFGENPATH=/nesi/nobackup/nesi00213/Pre-processing/SrfGen
export QCOREPATH=/nesi/nobackup/nesi00213/tmp/qcore
export GMPEMODELSPATH=/nesi/nobackup/nesi00213/tmp/post-processing/im_processing/computations/GMPE_models

# leave these alone
export PYTHONPATH=$QCOREPATH:$PYTHONPATH:$SRFGENPATH:$GMPEMODELSPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nesi/nobackup/nesi00213/gmt/gmt-stable/lib64
# to fix HDF5 file locking issue https://github.com/ALPSCore/ALPSCore/issues/348
export HDF5_USE_FILE_LOCKING=FALSE

echo "SRF generation starts"
date
python2 $SRFGENPATH/gcmt2srf.py $GCMT_FILE $SRFGENPATH/lp_generic1d-gp01_v1.vmod -n $NUMPROCS
echo "SRF generation completed"
date
echo "VM generation starts"
srun python2 $SRFGENPATH/srfinfo2vm.py "autosrf/*/Srf/*.info" -n $NUMPROCS
date
echo "VM generation completed"

## Only one VM is generated 3631755 as the other one (2122842) is in the ocean 

#!/bin/bash
# script version: slurm

# Please modify this file as needed, this is just a sample
#SBATCH --job-name=nhm2srf_vm
#SBATCH --account=nesi00213
#SBATCH --partition=NeSI
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --output nhm2srf_vm-%j.out
#SBATCH --error nhm2srf_vm-%j.err
###SBATCH --mail-type=all
###SBATCH --mail-user=test@test.com
###SBATCH --mem-per-cpu=16G
###SBATCH -C avx
#OpenMP+Hyperthreading works well for VM
###SBATCH --hint=nomultithread

## END HEADER
source machine_env.sh

# user specific parameters
export SRF_NUMPROCS=80
export VM_NUMPROCS=5
export OMP_NUM_THREADS=16

export SELECTION_FILE=`pwd`/nhm_selection_geo_si.txt #nhm_selection_sample2.txt #
# important paths
export SRFGENPATH=/nesi/nobackup/nesi00213/tmp/Pre-processing/SrfGen
export QCOREPATH=/nesi/nobackup/nesi00213/tmp/qcore
export GMPEMODELSPATH=/nesi/nobackup/nesi00213/tmp/post-processing/im_processing/computations/GMPE_models
export PATH=$PATH:/nesi/nobackup/nesi00213/tmp/Velocity-Model

# leave these alone
export PYTHONPATH=$QCOREPATH:$PYTHONPATH:$SRFGENPATH:$GMPEMODELSPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nesi/nobackup/nesi00213/gmt/gmt-stable/lib64

 #to fix HDF5 file locking issue https://github.com/ALPSCore/ALPSCore/issues/348
export HDF5_USE_FILE_LOCKING=FALSE

echo "SRF generation starts"
date
python2 $SRFGENPATH/NHM/nhm2srf.py $SELECTION_FILE -n $SRF_NUMPROCS
echo "SRF generation completed"
date
echo "VM generation starts"
srun python2 $SRFGENPATH/srfinfo2vm.py "autosrf/*/Srf/*.info" -n $VM_NUMPROCS
date
echo "VM generation completed"


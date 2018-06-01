#module swap Anaconda2/4.2.0-CrayGNU-2017.06 Anaconda2/5.0.1-CrayGNU-2017.06
module rm Anaconda2/4.2.0-CrayGNU-2017.06
module rm Anaconda2/5.0.1-CrayGNU-2017.06
source /nesi/transit/nesi00213/share/bashrc.uceq
#module rm Anaconda2/4.2.0-CrayGNU-2017.06
#module rm Anaconda2/5.0.1-CrayGNU-2017.06

# Replace with the actual Python module
module use /nesi/nobackup/nesi00213/spack/share/spack/modules/linux-sles12-x86_64/
module rm gcc/4.9.3
module load gcc/7.1.0
module load GDAL/2.2.1-CrayGNU-2017.06 
module load FFTW/3.3.7-CrayGNU-2017.06 
module load JasPer/2.0.12-CrayGNU-2017.06
module load curl-7.59.0-gcc-7.1.0-vkbeasr
module rm cray-hdf5/1.10.1.1
module rm HDF/4.2.13-CrayGNU-2017.06
module swap Anaconda2/4.2.0-CrayGNU-2017.06 Anaconda2/5.0.1-CrayGNU-2017.06
#module add cray-hdf5-parallel/1.10.1.1
export BINPROCESS=/nesi/transit/nesi00213/workflow/scripts
export PATH=$PATH:/nesi/nobackup/nesi00213/gmt/gmt-stable/bin:/nesi/transit/nesi00213/tools

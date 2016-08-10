# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = run_emod3d
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 6:00:00
#
# @ account_no = nesi00213
#
# @ output               = $(job_name).$(schedd_host).$(jobid).o
# @ error                = $(job_name).$(schedd_host).$(jobid).e
# @ notification         = never
# @ class                = General
# @ node = 16 
# @ tasks_per_node = 32 
# @ network.MPI = sn_all,not_shared,US
# @ task_affinity = core(1)
#
# @ queue
export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/20160809
source /nesi/projects/nesi00213/share/bashrc.uceq
python $BINPROCESS/set_runparams.py
date
poe /nesi/projects/nesi00213/EMOD3D/Mpi/Emod3d/V3.0.4/bin/powerpc-AIX-nesi2/emod3d-mpi -args "par=e3d.par"
date

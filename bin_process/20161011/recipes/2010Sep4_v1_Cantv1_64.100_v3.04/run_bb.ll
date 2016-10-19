# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = run_bb 
# @ account_no = nesi00213
# @ notification         = never
# @ class                = General
#

# @ step_name = hf
# @ job_type = serial 
# @ wall_clock_limit     = 1:00:00
# @ output = $(job_name).$(step_name).$(schedd_host).$(jobid).o
# @ error = $(job_name).$(step_name).$(schedd_host).$(jobid).e
# @ queue

# @ step_name = bb
# @ dependency = hf==0
# @ job_type = serial 
# @ wall_clock_limit     = 0:10:00
# @ output = $(job_name).$(step_name).$(schedd_host).$(jobid).o
# @ error = $(job_name).$(step_name).$(schedd_host).$(jobid).e
# @ queue
source ./version
export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/$VERSION

case $LOADL_STEP_NAME in
        hf)
                python $BINPROCESS/hfsims-stats.py ;;
        bb)
                python $BINPROCESS/match_seismo.py ;;
esac



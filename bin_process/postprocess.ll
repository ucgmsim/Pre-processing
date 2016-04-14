# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = postprocess 
#
#
# @ account_no = nesi00213
#
# @ output               = $(job_name).$(schedd_host).$(jobid).o
# @ error                = $(job_name).$(schedd_host).$(jobid).e
# @ notification         = never
# @ class                = General
#

# @ step_name = winbin_aio
# @ job_type = serial 
# @ wall_clock_limit     = 0:20:00
# @ queue

# @ step_name = gen_ts
# @ dependency = winbin_aio==0
# @ job_type = serial 
# @ wall_clock_limit     = 0:20:00
# @ queue

case $LOADL_STEP_NAME in
	winbin_aio)
		python winbin-aio.py ;;
	gen_ts)
		python gen_ts.py ;;
esac




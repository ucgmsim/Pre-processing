# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = postprocess 
# @ account_no = nesi00213
# @ notification         = never
# @ class                = General
#

# @ step_name = merge_tsP3
# @ job_type = parallel
# @ network.MPI = sn_all,shared,US
# @ task_affinity = core(1)
# @ node = 1
# @ tasks_per_node = 4
# @ wall_clock_limit    = 0:30:00
# @ output = $(job_name).$(step_name).$(schedd_host).$(jobid).o
# @ error = $(job_name).$(step_name).$(schedd_host).$(jobid).e
# @ queue
 
# @ step_name = winbin_aio
# @ job_type = serial 
# @ wall_clock_limit     = 0:20:00
# @ output = $(job_name).$(step_name).$(schedd_host).$(jobid).o
# @ error = $(job_name).$(step_name).$(schedd_host).$(jobid).e
# @ queue

# @ step_name = gen_ts
# @ dependency = merge_tsP3==0
# @ job_type = serial 
# @ wall_clock_limit     = 0:20:00
# @ output = $(job_name).$(step_name).$(schedd_host).$(jobid).o
# @ error = $(job_name).$(step_name).$(schedd_host).$(jobid).e
# @ queue

source ./version
export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/$VERSION

case $LOADL_STEP_NAME in
	merge_tsP3)
		exe=/nesi/projects/nesi00213/EMOD3D/merge_ts/merge_tsP3_par 
		ls -X LF/OutBin/*xyts-?????.e3d >tmp.filelist 
		NFILES=`cat tmp.filelist |wc -l|sed -e 's/ //g'` 
		OUTFILE=`head -1 tmp.filelist  |sed -e 's/-[0-9]*.e3d/.e3d/g'` 
		echo "NFILES=$NFILES OUTFILE=$OUTFILE" 
		time poe $exe filelist=tmp.filelist outfile=$OUTFILE nfiles=$NFILES ;;
	winbin_aio)
		python $BINPROCESS/winbin-aio.py ;;
	gen_ts)
		python $BINPROCESS/gen_ts.py ;;
esac




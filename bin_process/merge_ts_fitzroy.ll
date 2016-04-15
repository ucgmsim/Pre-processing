# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = merge 
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 2:00:00
#
# @ account_no = nesi00213
#
# @ network.MPI = sn_all,shared,US
# @ task_affinity = core(1)
#
# @ output               = $(job_name).$(schedd_host).$(jobid).o
# @ error                = $(job_name).$(schedd_host).$(jobid).e
# @ notification         = never
# @ class                = General
# @ node = 1
# @ tasks_per_node = 4 
#
# @ queue

home=/gpfs_external/filesets/nesi/home/pletzera
exe=/nesi/projects/nesi00213/EMOD3D/merge_ts/merge_tsP3_par

ls -1 OutBin/*xyts-?????.e3d >tmp.filelist
export NFILES=`cat tmp.filelist |wc -l`
export OUTFILE=`head -1 tmp.filelist  |sed -e 's/-[0-9]*.e3d/.e3d/g'`

echo "NFILES=$NFILES OUTFILE=$OUTFILE"

time poe $exe filelist=tmp.filelist outfile=$OUTFILE nfiles=$NFILES




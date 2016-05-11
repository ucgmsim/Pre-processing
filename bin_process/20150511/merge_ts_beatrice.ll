# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = merge
#
# @ job_type = MPICH 
# @ node = 1 
# @ tasks_per_node = 4 
#
# @ wall_clock_limit = 04:00:00
#
# Groups to select from: UC, UC_merit, NZ, NZ_merit
# @ group = NZ_merit 
# Your project number, either bfcs or nesi, followed by 5 digits
# @ account_no = nesi00213
# @ output = $(job_name).$(schedd_host).$(jobid).out
# @ error = $(job_name).$(schedd_host).$(jobid).err
# To receive an email when yourjob has completed:
# @ notification       = never
# @ notify_user        = myemail@gmail.com
 
# @ class = p7linux
#
# To improve performance it is important to specify task and memory affinities
# @ rset = rset_mcm_affinity
## @ task_affinity = core(1) 
# @ network.MPI_LAPI = sn_single,shared,US,,instances=2
# @ environment = COPY_ALL
# @ queue
  
# suggested environment settings:
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM
 
# Launch parallel job
module load openmpi/1.10_adv
module load gcc/adv70
ls -X OutBin/*xyts-?????.e3d >tmp.filelist
export NFILES=`cat tmp.filelist |wc -l|sed -e 's/ //g'`
export OUTFILE=`head -1 tmp.filelist  |sed -e 's/-[0-9]*.e3d/.e3d/g'`_beatrice
echo $NFILES $OUTFILE
mpirun /hpc/scratch/nesi00213/EMOD3D/merge_ts/bin/ppc64-Linux-p2n14-c-gcc/merge_tsP3_par filelist=tmp.filelist outfile=$OUTFILE nfiles=$NFILES


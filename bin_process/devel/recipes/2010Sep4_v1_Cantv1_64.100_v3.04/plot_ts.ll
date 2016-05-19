# Example MPI LoadLeveler Job file
# @ shell = /bin/bash
#
# @ job_name = plot_ts
#
# @ job_type = serial 
#
# @ wall_clock_limit = 00:30:00
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
 
# @ class = p7linux_serial
#
# To improve performance it is important to specify task and memory affinities
# @ rset = rset_mcm_affinity
# @ environment = COPY_ALL
# @ queue
  
# suggested environment settings:
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM
export BINPROCESS=/hpc/scratch/nesi00213/Pre-processing/bin_process/devel
bash $BINPROCESS/plot_ts.sh
bash $BINPROCESS/make_movie.sh

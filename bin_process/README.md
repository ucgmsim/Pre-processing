How to run a emod3d workflow on NIWA Fitzroy

Log into Fitzroy
```
ssh <username>@fitzroy.nesi.org.nz
```
If you don't have an account, contact hpcf-admins@niwa.co.nz .

Go to the project directory
```
cd /nesi/projects/nesi00213/RunFolder
```

Make a directory with a descriptive name
```
mkdir LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04
```

Go into the directory
```
cd LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04
```

This repository contains a number of recipes that can be reused. See if a recipe for your modelling is already available

```
ls /nesi/projects/nesi00213/Pre-processing/bin_process/stable/recipes
```

If you found a recipe, copy them to this directory. 
If there is no recipe available, select one that looks most relevant to the model you will be running, and edit them to suit your needs.

```
cp /nesi/projects/nesi00213/Pre-processing/bin_process/stable/recipes/2010Sep4_v1_Cantv1_64.100_v3.04/* .
```

Examine/Edit params.py 

Examine/Edit run_emod3d.ll. Pay attention to the wall_clock_limit, node, and tasks_per_node. 
Each node can run 32 tasks in normal mode. You can also specify 64 tasks per node in SMT mode.

```
# @ shell = /bin/bash
#
# @ job_name = run_emod3d
#
# @ job_type = parallel
#
# @ wall_clock_limit     = 4:00:00
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
export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/20160511
python $BINPROCESS/set_runparams.py
poe /nesi/projects/nesi00213/EMOD3D/Mpi/Emod3d/V3.0.4/bin/powerpc-AIX-nesi2/emod3d-mpi -args "par=e3d.par"

```

Submit the job.

```
llsubmit run_emod3d.ll
```

Examine/edit post_emod3d.ll. This has a number of steps organized to execute in series/parallel. 

(1). merge_tsP3 -> gen_ts 
(2). winbin_aio -> hf -> bb

Except for merge_tsP3, all steps are single process jobs. Just pay attention to wall_clock_limit. 
The specified wall_clock_limit's have been tested sufficient.

```
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

# @ step_name = hf
# @ dependency = winbin_aio==0
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

export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/20160511

case $LOADL_STEP_NAME in
	merge_tsP3)
		exe=/nesi/projects/nesi00213/EMOD3D/merge_ts/merge_tsP3_par 
		ls -X OutBin/*xyts-?????.e3d >tmp.filelist 
		NFILES=`cat tmp.filelist |wc -l|sed -e 's/ //g'` 
		OUTFILE=`head -1 tmp.filelist  |sed -e 's/-[0-9]*.e3d/.e3d/g'` 
		echo "NFILES=$NFILES OUTFILE=$OUTFILE" 
		time poe $exe filelist=tmp.filelist outfile=$OUTFILE nfiles=$NFILES ;;
	winbin_aio)
		python $BINPROCESS/winbin-aio.py ;;
	gen_ts)
		python $BINPROCESS/gen_ts.py ;;

	hf)
		python $BINPROCESS/hfsims-stats.py ;;
	bb)
		python $BINPROCESS/match_seismo.py ;;
esac
```

Submit the post-process job
```
llsubmit post_emod3d.ll
```
After the job is completed, you should see a merged _xtys.e3d file under "OutBin", newly created "Vel" and "TSlice" under current directory.
Outputs for hf and bb modelling are stored in newly created (if previously non-existent) HF_xxxxx and BB_xxxxx (where xxxxx are replaced by hf_run_name and bb_run_name specified in params.py) under RunFolder.

plot_ts needs to be done by a POWER-based linux. Go to the linux node, and edit plot_ts.ll if required, then submit

```
llsubmit plot_ts.ll
```

If a movie file from generated Png files is desired, go to TSlice/Png directory and run

```
convert -delay 5 *.png <filename>.gif
```


If everything went smoothly, contribute to the recipe repository!

Known issues:
1. If a step in post_emod3d.ll has an error and terminates, and a job needs to be resubmitted, steps that had no problem will be computed again.
To avoid this, you will need to edit post_emod3d.ll before resubmission.
2. HF produces Acc result and you can use acc2vel.py to convert to Vel (groundMotionStationAnalysis needs to Vel input). 
As acc2vel.py uses parameters set in params.py, you will need to run it from LP.... directory or need to create a symbolic link to params.py from the HF directory and run it there.


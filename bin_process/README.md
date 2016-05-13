# How to run a emod3d workflow on NIWA Fitzroy

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

Examine/Edit params.py. 
This file is roughly composed of 7 parts. Some important ones that need more attentions are listed below.
```
........
######## Global CONSTANTS ########
run_name = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04'
version = '3.0.4'
# LPSIM directory name under RunFolder
extended_run_name = 'LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04'
# HFSIM directory name under RunFolder
hf_run_name = 'HFSim-2010Sept4b_Cant1D_v2-v5.4.4-rvf0.8_dt'
bb_run_name = 'BBSim-2010Sept4b_Cantv1_64-h0.100_dt_Vs30_500'
# prefix used for naming OutBin/prefix{_xyts.e3d,_seis.e3d}
output_prefix = run_name #If output file names are different from the directory name. 
....
# keep as int, processed before writing to file
# only product is stored
n_proc_x = 8
n_proc_y = 8
n_proc_z = 8

### folowing values are dependent on the velocity model used
# low pass frequency, filter for output seimograms
flo = '1.0' # hertz
# spatial grid spacing
hh = '0.100' # km
# x, y, z grid size (multiples of grid spacing)
nx = '1400'
ny = '1200'
nz = '460'
# model reference location
MODEL_LAT = '-43.6000'
MODEL_LON = '172.3000'
MODEL_ROT = '-10.0'
###

# cap number of timesteps in simulation, not all timesteps have outputs
# max simulation timeperiod = nt * dt eg: 10,000 * 0.005 = 50 seconds
nt = '20000'
# dt should be 0.005 (or smaller), small increments required in simulation
dt = '0.005'
# how often to save outputs (measured in simulation timesteps)
DUMP_ITINC = '4000' # nt

....

# which time slices to iterate over
ts_start = '0'     # first one is 0
ts_inc = '1'       # increment, larger than 1 to skip
ts_total = '400'   # number of slices to generate. sim time = ts_total * dt * dt_ts
....
vel_mod_params_dir = os.path.join(user_root, 'VelocityModel/ModelParams')
vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel', v_mod_ver)
....
# files
srf_file = os.path.join(srf_dir, '2010Sept4_m7pt1/Srf/bev01.srf')
stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_nz01-h0.100.statcords')


############# winbin-aio ##############
....
############### gen_ts ###################
....
################## plot_ts ####################
....
############ acc2vel and match_seismo #############
....
########### gen_statgrid ##########
....

```

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

This will produce a number of xxxx_xyts.e3d and xxxx_seis.e3d files under OutBin.

Examine/edit post_emod3d.ll. This has a number of steps organized to execute in series/parallel. 

+ merge_tsP3 -> gen_ts 
+ winbin_aio -> hf -> bb

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


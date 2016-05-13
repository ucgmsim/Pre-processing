# How to run QuakeCoRE workflow on NIWA Fitzroy

Log into Fitzroy. If you don't have an account, contact hpcf-admins@niwa.co.nz
```
ssh <username>@fitzroy.nesi.org.nz
```

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
This file is roughly composed of 7 parts. Some important ones that need more attentions are listed below. In particular, the product of n_proc_x, n_proc_y and n_proc_z should be the same as the total number of processes. (512 in this example)
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
Each node can run 32 tasks in normal mode. You can also specify 64 tasks per node in SMT mode. In this example, we use node * tasks_per_node = 512 tasks.

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

This will produce a number of xxxx_xyts-????.e3d and xxxx_seis-????.e3d files under OutBin. 

Examine/edit post_emod3d.ll. This has a number of steps organized to execute in series/parallel. 

+ merge_tsP3 -> gen_ts 
+ winbin_aio -> hf -> bb

Except for merge_tsP3, all steps are single process jobs (merge_tsP3 uses 4 processes which should be adequate for most cases). You may need to edit wall_clock_limit for each step (most likely unnecessary)
Note that merge_tsP3 assumes that all the xxxx_xyts-????.e3d files are kept in **OutBin** directory, and the merged file will be also stored there.
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
Output of each steps is:
+ merge_tsp3: A merged file LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_xyts.e3d under "OutBin"
+ winbin-aio: A newly created "Vel" directory and files that look like below. These are produced from OutBin/xxxx_seis-????.e3d files.
```
ADCS.000  CACS.ver  CHHC.090  CSHS.000  D06C.ver  D15C.090  DSLC.000  GDLC.ver  HHSS.090  HVSC.000  KPOC.ver  LPOC.090  MTHS.000  NNBS.ver  PPHS.090  RHSC.000  ROLC.ver  SHLC.090  SPFS.000  STKS.ver  WCSS.090
ADCS.090  CBGS.000  CHHC.ver  CSHS.090  D09C.000  D15C.ver  DSLC.090  GODS.000  HHSS.ver  HVSC.090  LINC.000  LPOC.ver  MTHS.090  ...
```
+ gen_ts: A newly created TSlice/TSFiles directory  and LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04_ts0000.0 ~ _ts0399.0 files. (399 is due to ts_total = 400 in params.py). This step depends on merge_tsP3.
+ hf : RunFolder/HFSim-2010Sept4b_Cant1D_v2-v5.4.4-rvf0.8_dt directory (if previously non-existent) and Acc subdirectory
+ bb:  RunFolder/BBSim-2010Sept4b_Cantv1_64-h0.100_dt_Vs30_500 directory (if previously non-existent) and Vel and Acc subdirectory 

Go to Linux node. Edit plot_ts.ll (Can be used untouched in most cases), and submit
```
llsubmit plot_ts.ll
```
This will produce Png and Ps files under TSlice directory per each time slice.
In addition, plot_ts.ll execute make_movie.sh after plotting is completed, which will produce a movie named after the run_name in the working directory.

```
-rw-r--r-- 1 seb56 n00213 52636372 May 13 15:38 LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04.gif
```

If everything went smoothly, contribute to the recipe repository!

Known issues:

1. If a step in post_emod3d.ll has an error and terminates, and a job needs to be resubmitted, steps that had no problem will be computed again. To avoid this, you will need to edit post_emod3d.ll before resubmission.
2. HF produces Acc result and you can use acc2vel.py to convert to Vel (groundMotionStationAnalysis needs to Vel input). 
As acc2vel.py uses parameters set in params.py, you will need to run it from LP.... directory or need to create a symbolic link to params.py from the HF directory and run it there.
3. HF Acc and LP Vel are used as inputs for BB, and historic garbage files (if exists) may impact the output. These directories are now **all emptied before the computation**. If files need to be kept, **make a backup** manually.
4. Currently, plot_ts.ll needs to be submitted at Beatrice and it requires editing of path in e3d.par.

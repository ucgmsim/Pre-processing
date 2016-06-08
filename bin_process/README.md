# How to run QuakeCoRE workflow on NIWA Fitzroy

Version. 20160607

## Quick Summary ##


1. Go to /nesi/projects/nesi00213/RunFolder
2. ./install.sh and answer questions
3. New directory will be created under your subdirectory
4. Edit params.py
5. llsubmit run_emod3d.ll
5. llsubmit post_emod3d.ll
6. For HF/BB: install_bb.sh
7. For HF/BB: llsubmit run_bb.ll
6. For plotting: (Linux) bash plot_and_ani.sh


## Step-by-step Instruction ##

Log into Fitzroy. If you don't have an account, contact hpcf-admins@niwa.co.nz
```
ssh <username>@fitzroy.nesi.org.nz
```

Go to the project directory
```
cd /nesi/projects/nesi00213/RunFolder
```

Run install.py
```
baes@nesi2/nesi/projects/nesi00213/RunFolder >./install.sh
****************************************************************************************************
                                     EMOD3D Job Preparationi Ver.20160607
****************************************************************************************************
====================================================================================================
Select Rupture Model - Step 1.
====================================================================================================
 1. 2010Dec26_m4pt7
 2. 2010Oct19_m4pt8
 3. 2010Sept4_m7pt1
 4. 2011Apr16_m5pt0
 5. 2011Apr16_m5pt0_old
 6. 2011Dec23a_m5pt8
 7. 2011Dec23a_m5pt8_old
 8. 2011Dec23b_m5pt9
 9. 2011Dec23b_m5pt9_old
10. 2011Feb22_m6pt2
11. 2011Feb22_m6pt2_
12. 2011Feb22_m6pt2_old
13. 2011June13a_m5pt3
14. 2011June13b_m6pt0
15. 2011June13b_m6pt0_old
16. 2011June21_m5pt2
17. 2011June21_m5pt2_old
18. 2016Feb14_m5pt8
19. Multi_Segment
20. Multi_Segment_v1
21. RupModelDevelopment
Enter the number you wish to select (1-21):3
2010Sept4_m7pt1
====================================================================================================
Select Rupture Model - Step 2.
====================================================================================================
 1. bev01.srf
 2. bev01_s103245.srf
 3. bev01_s103245v1.srf
 4. bev01_s103245v2.srf
 5. bev01_s103246.srf
 6. bev01_s103246v1.srf
 7. bev01_s103246v2.srf
 8. bev01_s103247.srf
 9. bev01_s103247v1.srf
10. bev01_s103247v2.srf
11. bev01_s103248.srf
12. bev01_s103248v1.srf
13. bev01_s103248v2.srf
14. bev01_s103249.srf
15. bev01_s103249v1.srf
16. bev01_s103249v2.srf
17. bev01_s103250.srf
18. bev01_s103250v1.srf
19. bev01_s103250v2.srf
20. bev01_s103251.srf
21. bev01_s103251v1.srf
22. bev01_s103251v2.srf
23. bev01_s103252.srf
24. bev01_s103252v1.srf
25. bev01_s103252v2.srf
26. bev01_s103253.srf
27. bev01_s103253v1.srf
28. bev01_s103253v2.srf
29. bev01_s103254.srf
30. bev01_s103254v1.srf
31. bev01_s103254v2.srf
Enter the number you wish to select (1-31):1
bev01.srf
Corresponding Stock file is also found:
/nesi/projects/nesi00213/RupModel/2010Sept4_m7pt1/Stoch/bev01.stoch
====================================================================================================
Select HH 
====================================================================================================
 1. 0.100
 2. 0.200
 3. 0.400
Enter the number you wish to select (1-3):3
0.400
====================================================================================================
Select one of available CanterburyVelocityModels (from /nesi/projects/nesi00213/CanterburyVelocityModel)
====================================================================================================
 1. v1.01remote
 2. v1.02
 3. v1.64
Enter the number you wish to select (1-3):3
v1.64
/nesi/projects/nesi00213/CanterburyVelocityModel/v1.64
====================================================================================================
Automated Name:  2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606
====================================================================================================
Do you wish to proceed?
 1. Yes
 2. No
Enter the number you wish to select (1-2):1
Yes?  True
====================================================================================================
Choose one of available recipes (from /nesi/projects/nesi00213/Pre-processing/bin_process/20160607/recipes)
====================================================================================================
 1. 2010Sep4_v1_Cantv1_64.100_v3.04
 2. 2010Dec26_m4.7pt_v1_Cant1D_v1-nz01-h0.500_V3.04
Enter the number you wish to select (1-2):1
2010Sep4_v1_Cantv1_64.100_v3.04
****************************************************************************************************
To be created: 
2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606
Recipe to be copied from 
/nesi/projects/nesi00213/Pre-processing/bin_process/20160607/recipes/2010Sep4_v1_Cantv1_64.100_v3.04
****************************************************************************************************
Do you wish to proceed?
 1. Yes
 2. No
Enter the number you wish to select (1-2):1
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606 : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/LF : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/HF : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/BB : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/Figures : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/run_bb.ll : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/run_emod3d.ll : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/params.py : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/post_emod3d.ll : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/params_base.py : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/install_bb.py : 750
Installation completed
====================================================================================================
Instructions
====================================================================================================
    1.   cd /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606
    2.   Edit params.py
    3.   llsubmit run_emod3d.ll
    4.   llsubmit post_emod3d.ll
    5.   (Linux) bash plot_ts.sh
    6.   python install_bb.py
    7.   llsubmit run_bb.ll

```

Examine/Edit params.py. 
This file is roughly composed of 6 parts. Some important ones that need more attentions are listed below. In particular, the product of n_proc_x, n_proc_y and n_proc_z should be the same as the total number of processes. (512 in this example)
```
........
######## Global CONSTANTS ########
........
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
export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/20160607
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


export BINPROCESS=/nesi/projects/nesi00213/Pre-processing/bin_process/20160607

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

esac
```

Submit the post-process job
```
llsubmit post_emod3d.ll
```
Output of each steps is:
+ merge_tsp3: A merged file ???_xyts.e3d under "OutBin"
+ winbin-aio: A newly created "Vel" directory and files that look like below. These are produced from OutBin/xxxx_seis-????.e3d files.
```
ADCS.000  CACS.ver  CHHC.090  CSHS.000  D06C.ver  D15C.090  DSLC.000  GDLC.ver  HHSS.090  HVSC.000  KPOC.ver  LPOC.090  MTHS.000  NNBS.ver  PPHS.090  RHSC.000  ROLC.ver  SHLC.090  SPFS.000  STKS.ver  WCSS.090
ADCS.090  CBGS.000  CHHC.ver  CSHS.090  D09C.000  D15C.ver  DSLC.090  GODS.000  HHSS.ver  HVSC.090  LINC.000  LPOC.ver  MTHS.090  ...
```
+ gen_ts: A newly created TSlice/TSFiles directory  and 2010Sept4_v1_Cantv1_64-h0.100_v3.04_ts0000.0 ~ _ts0399.0 files. (399 is due to ts_total = 400 in params.py). This step depends on merge_tsP3.
+ 

Run install_bb.py to prepare for HF/BB runs
```
baes@nesi2/nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606 >install_bb.sh

****************************************************************************************************
                                     EMOD3D HF/BB Preparationi Ver.20160607
****************************************************************************************************
====================================================================================================
Select one of 1D Velocity models (from /nesi/projects/nesi00213/VelocityModel/Mod-1D)
====================================================================================================
 1. /nesi/projects/nesi00213/VelocityModel/Mod-1D/Cant1D_v1-midQ.1d
 2. /nesi/projects/nesi00213/VelocityModel/Mod-1D/Cant1D_v1.1d
 3. /nesi/projects/nesi00213/VelocityModel/Mod-1D/Cant1D_v2-midQ.1d
Enter the number you wish to select (1-3):1
/nesi/projects/nesi00213/VelocityModel/Mod-1D/Cant1D_v1-midQ.1d
====================================================================================================
- Vel. Model 1D: Cant1D_v1-midQ
- hf_sim_bin: hb_high_v5.4.5
- hf_rvfac: 0.8
- hf_sdrop: 50
- hf_kappa: 0.045
====================================================================================================
Automated Name:  Cant1D_v1-midQ_hfv5p4p5_rvf0p8_sd50_k0p045
====================================================================================================
Do you wish to proceed?
 1. Yes
 2. No
Enter the number you wish to select (1-2):1
Yes?  True
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/HF/Cant1D_v1-midQ_hfv5p4p5_rvf0p8_sd50_k0p045 : 750
Permission /nesi/projects/nesi00213/RunFolder/baes/2010Sept4_bev01_VMv1p64-h0p400_EMODv3p0p4_160606/BB/Cant1D_v1-midQ_hfv5p4p5_rvf0p8_sd50_k0p045 : 750
```

and submit the job (after editing run_bb.ll if necessary)

```
llsubmit run_bb.ll
```


Go to Linux node. Edit plot_ts.ll (Can be used untouched in most cases).
```
plot_and_ani.sh
```
This will produce Png and Ps files under TSlice directory per each time slice.
This also makes a movie file in the current directory.

```
-rw-r--r-- 1 seb56 n00213 52636372 May 13 15:38 LPSim-2010Sept4_v1_Cantv1_64-h0.100_v3.04.gif
```

If everything went smoothly, contribute to the recipe repository!

Known issues:

1. If a step in post_emod3d.ll has an error and terminates, and a job needs to be resubmitted, steps that had no problem will be computed again. To avoid this, you will need to edit post_emod3d.ll before resubmission.
2. HF Acc and LP Vel are used as inputs for BB, and historic garbage files (if exists) may impact the output. These directories are now **all emptied before the computation**. If files need to be kept, **make a backup** manually.


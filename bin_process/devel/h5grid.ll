# extracts velocities from seis files for each grid point
# generates HF velocities for each grid point
# creates BB data
# stores BB in HDF5

# PRIOR TO RUNNING THIS
# gen_statgrid.py - create station points
# lf sim job - with additional station file above
# RUNNING THIS WILL:
# A: Extract seis data from LF output at virtual stations to H5 file
# B: Run HF sim for virtual stations.
# C: Run seis match for h5 LF, save next to HF.
# D: store processed data to H5 (serial job)

# @shell = /bin/bash
# @job_name = h5_grid
# @output = $(job_name).out
# @error = $(job_name).err
# @account_no = nesi00213
# @class = General
# @notification = never

# LF SEIS -> HDF5
# @step_name = lf_extract
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# @wall_clock_limit = 0:30:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

# HF SIM -> Binary Data Files
# nt=20000, 1400x1200 takes ~30 minutes
# @step_name = hf_sim
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# @wall_clock_limit = 1:00:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

# SEIS MATCH -> Binary intermediate
# @step_name = bb_match
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# @dependency = (hf_sim==0) && (lf_extract==0)
# probably doesn't need to be this long, this is roughly single threaded laptop time
# @wall_clock_limit = 1:00:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

# BB WRITE -> HDF5
# @step_name = bb_write
# @dependency = (bb_match==0)
# @wall_clock_limit = 0:30:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

source /nesi/projects/nesi00213/share/bashrc.uceq
export PYTHONPATH=$PYTHONPATH:$(pwd)

case $LOADL_STEP_NAME in
    lf_extract)
        echo "lf extract starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./h5lf.py
        ;;
    hf_sim)
        echo "hf_sim starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./hfsims-stats.py 64 EMOD3D/StochSim/Src/V5.4/hb_high_v5.4.5_binmod
        ;;
    bb_match)
        echo "bb_match starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./h5match.py 64
        ;;
    bb_write)
        echo "bb_write starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./h5bb.py
        ;;
esac

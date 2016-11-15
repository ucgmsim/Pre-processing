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
# @job_name = seisfind_store
# @output = $(job_name).out
# @error = $(job_name).err
# @account_no = nesi00213
# @class = General
# @notification = never

# HF SIM -> Binary Data Files
# nt=20000, 1400x1200 takes ~30 minutes
# @step_name = hf_sim
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# @wall_clock_limit = 4:00:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

# LF SEIS -> HDF5
# @step_name = lf2bb
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# @wall_clock_limit = 10:00:00
# @output = $(job_name).$(step_name).o
# @error = $(job_name).$(step_name).e
# @queue

source /nesi/projects/nesi00213/share/bashrc.uceq
export PYTHONPATH=$PYTHONPATH:$(pwd)

case $LOADL_STEP_NAME in
    hf_sim)
        date
        echo "hf_sim starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./hfsims-stats.py 64 EMOD3D/StochSim/Src/V5.4/hb_high_v5.4.5_binmod
        date
        echo "hf_sim finished."
        ;;
    lf2bb)
        date
        echo "lf2bb starting..."
        /opt/niwa/Python/AIX/2.7.5/bin/python ./seisfind_store.py 64
        date
        echo "lf2bb finished."
        ;;
esac

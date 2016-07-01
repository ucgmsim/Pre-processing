# extracts velocities from seis files for each grid point
# generates HF velocities for each grid point
# creates BB data
# stores BB in HDF5

# @shell = /bin/bash
# @job_name = hf_grid
# @job_type = parallel
# @node = 1
# @tasks_per_node = 64
# nt=20000, 1400x1200 takes ~30 minutes
# @wall_clock_limit = 0:45:00
# @account_no = nesi00213
# @output = output
# @error = error
# @notification = never
# @class = General
# @queue

# only HF generation so far
cd /nesi/projects/nesi00213/scratch/vap30
/opt/niwa/Python/AIX/2.7.5/bin/python hfsims-stats.py 64

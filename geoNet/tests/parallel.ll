# Example MPI LoadLeveler Job file
# @shell = /bin/bash
#
# @job_name = process_data
# @account_no = nesi00213
# @notification         = never
# @class                = General
#

# @job_type = serial
# @wall_clock_limit     = 3:00:00
# @output = $(job_name).$(schedd_host).$(jobid).o
# @error = $(job_name).$(schedd_host).$(jobid).e
# @queue

source ~/.bashrc
#export ANALYSIS="$RUN"/groundMotionStationAnalysis
/opt/niwa/Python/AIX/2.7.5/bin/python Mw4pt1_20110707_processData.py

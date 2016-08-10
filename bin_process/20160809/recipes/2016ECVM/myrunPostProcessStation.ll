# @shell = /bin/bash
# @job_name = computeIM
# @job_type = parallel
# @node = 1
# @tasks_per_node = 30
# @wall_clock_limit = 0:10:00
# @account_no = nesi00213
# @output = $(job_name).out
# @error = $(job_name).err
# @notification = never
# @class = General
# @queue
source /home/ahsan/.bashrc
#export ANALYSIS=$RUN/groundMotionStationAnalysis
export ANALYSIS=$QC/groundMotionStationAnalysis

# Launch parallel job
tasks_per_node=30   #this should automatically be read from above somehow!!!


#/opt/niwa/Python/AIX/2.7.5/bin/python $ANALYSIS/runPostProcessStation.py parametersStation.py $tasks_per_node
#parametersStation_Mw4pt6Sep2010.py

/opt/niwa/Python/AIX/2.7.5/bin/python $ANALYSIS/runPostProcessStation.py parametersStation_Mw4pt6Nov2010 $tasks_per_node

 

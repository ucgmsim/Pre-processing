bin_process_path='/nesi/projects/nesi00213/Pre-processing/bin_process'
import glob
import os.path
import sys
from version import *
bin_process_dir = os.path.join(bin_process_path,bin_process_ver)
sys.path.append(os.path.abspath(os.path.curdir))
from params import *
from params_bb_base import *

sys.path.append(bin_process_dir)
from shared import *

def confirm(q):
    show_horizontal_line
    print q
    return show_yes_no_question()


hf_sim_dirs=glob.glob('HF/%s/*'%hf_sim_basedir)

f_template=open('run_bb.ll.template')
template=f_template.readlines()
template
str_template=''.join(template)

submit_yes = confirm("Also submit the job for you?")

for hf_sim_dir in hf_sim_dirs:
    bb_sim_dir = hf_sim_dir.replace('HF','BB')
    
    str=str_template.replace("$hf_sim_dir",hf_sim_dir)
    str=str_template.replace("$bb_sim_dir",bb_sim_dir)

    variation = '_'.join(hf_sim_dir.split('/')[1:3])
    fname_llscript='run_bb_%s.ll'%variation
    f_llscript=open(fname_llscript,'w')
    f_llscript.write('# script version: %s\n'%bin_process_ver)
    f_llscript.write(str)
    f_llscript.close()
    print "Loadleveler script %s written" %fname_llscript
    if submit_yes:
        print "Submitting %s" %fname_llscript
        res=execute_cmd("llsubmit %s"%fname_llscript)
        print res
    else:
        print "User chose to submit the job manually"


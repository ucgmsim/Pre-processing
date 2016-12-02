bin_process_path='/nesi/projects/nesi00213/Pre-processing/bin_process'
import glob
import os.path
import sys
from version import *
bin_process_dir = os.path.join(bin_process_path,bin_process_ver)
sys.path.append(bin_process_dir)
from shared import *

def confirm(q):
    show_horizontal_line
    print q
    return show_yes_no_question()

glob.glob('LF/*')
lf_sim_dirs=glob.glob('LF/*')
f_template=open('post_emod3d.ll.template')
template=f_template.readlines()
template
str_template=''.join(template)

submit_yes = confirm("Also submit the job for you?")

for lf_sim_dir in lf_sim_dirs:
    str=str_template.replace("$lf_sim_dir",lf_sim_dir)
    rup_mod = lf_sim_dir.split('/')[1]
    fname_llscript='post_emod3d_%s.ll'%rup_mod
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


bin_process_path='/nesi/projects/nesi00213/Pre-processing/bin_process'
import glob
import os.path
import sys
from version import *
bin_process_dir = os.path.join(bin_process_path,bin_process_ver)
sys.path.append(bin_process_dir)
import shared

execfile(os.path.join(bin_process_dir,"set_runparams.py"))
glob.glob('LF/*')
lf_sim_dirs=glob.glob('LF/*')
f_template=open('run_emod3d.ll.template')
template=f_template.readlines()
template
str_template=''.join(template)
for lf_sim_dir in lf_sim_dirs:
    str=str_template.replace("$lf_sim_dir",lf_sim_dir)
    rup_mod = lf_sim_dir.split('/')[1]
    fname_llscript='run_emod3d_%s.ll'%rup_mod
    f_llscript=open(fname_llscript,'w')
    f_llscript.write('# script version: %s\n'%bin_process_ver)
    f_llscript.write(str)
    f_llscript.close()
    res=shared.execute_cmd("llsubmit %s"%fname_llscript)
    print res

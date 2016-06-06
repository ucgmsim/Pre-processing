import os.path
import sys
sys.path.append(os.path.abspath(os.path.curdir))
from params import *
from params_base import *
import glob
from shared import *

bin_process_ver = 'devel'
  
def q1(v_mod_1d_dir): 
    show_horizontal_line()
    print "Select one of 1D Velocity models (from %s)" %v_mod_1d_dir
    show_horizontal_line()
 
    v_mod_1d_options = glob.glob(os.path.join(v_mod_1d_dir,'*.1d'))
    v_mod_1d_options.sort()

    v_mod_1d_selected = show_multiple_choice(v_mod_1d_options)
    print v_mod_1d_selected #full path
    v_mod_1d_name = os.path.basename(v_mod_1d_selected).replace('.1d','')
#    print v_mod_1d_name

    return v_mod_1d_name,v_mod_1d_selected

def q2(v_mod_1d_name):
    hfVString='hf'+os.path.basename(hf_sim_bin).split('_')[-1]
    hf_run_name=v_mod_1d_name+'_'+hfVString+'_rvf'+hf_rvfac+'_sd'+hf_sdrop+'_k'+hf_kappa
    hf_run_name=hf_run_name.replace('.','p')
    show_horizontal_line()
    print "- Vel. Model 1D: %s" %v_mod_1d_name
    print "- hf_sim_bin: %s" %os.path.basename(hf_sim_bin)
    print "- hf_rvfac: %s" %hf_rvfac
    print "- hf_sdrop: %s" %hf_sdrop
    print "- hf_kappa: %s" %hf_kappa
    yes = confirm_name(hf_run_name)
    return yes, hf_run_name
    
def action(hf_run_name,v_mod_1d_selected):
    hf_sim_dir, bb_sim_dir = os.path.join(hf_dir,hf_run_name), os.path.join(bb_dir,hf_run_name)
    dir_list =  [hf_sim_dir, bb_sim_dir]
    
    verify_user_dirs(dir_list)
    
    f=open(os.path.join(sim_dir,"params_base_bb.py"),"w");
    f.write("hf_run_name='%s'\n"%hf_run_name)
    f.write("hf_sim_dir='%s'\n"%hf_sim_dir)
    f.write("bb_sim_dir='%s'\n"%bb_sim_dir)
    f.write("hf_accdir='%s'\n"%os.path.join(hf_sim_dir,"Acc"))
    f.write("hf_veldir='%s'\n"%os.path.join(hf_sim_dir,"Vel"))
    f.write("bb_accdir='%s'\n"%os.path.join(bb_sim_dir,"Acc"))
    f.write("bb_veldir='%s'\n"%os.path.join(bb_sim_dir,"Vel"))
    f.write("hf_v_model='%s'\n"%v_mod_1d_selected)
    f.close()

    #creating symbolic link between matching HF and BB directories
#    try: 
#        os.symlink(bb_sim_dir,os.path.join(hf_sim_dir,"BB"))
#    except OSError:
#        print "Directory already exists: %s" %os.path.join(hf_sim_dir,"BB")
#    try: 
#        os.symlink(hf_sim_dir,os.path.join(bb_sim_dir,"HF"))
#    except OSError:
#        print "Directory already exists: %s" %os.path.join(bb_sim_dir,"HF")

    set_permission(hf_sim_dir)
    set_permission(bb_sim_dir)


def main():
    show_horizontal_line(c="*")
    print " "*37+"EMOD3D HF/BB Preparationi Ver."+bin_process_ver
    show_horizontal_line(c="*")

    v_mod_1d_name, v_mod_1d_selected = q1(v_mod_1d_dir)
    yes, hf_run_name = q2(v_mod_1d_name)
    hf_run_name = add_name_suffix(hf_run_name,yes)
    action(hf_run_name,v_mod_1d_selected)
    
    

if __name__ == "__main__":
    main() 

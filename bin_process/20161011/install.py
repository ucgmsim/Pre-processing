#!/usr/bin/env python2
emod3d_version = '3.0.4'

global_root='/nesi/projects/nesi00213/'
bin_process_path = '/nesi/projects/nesi00213/Pre-processing/bin_process/'
install_bb_name='install_bb.py'

params_vel = 'params_vel.py'
params_src = 'params_src.py'

import os
import os.path
import sys
import datetime
import glob
import shutil
import getpass
import importlib

from version import *

bin_process_dir = os.path.join(bin_process_path,bin_process_ver)
sys.path.append(bin_process_dir)  #unnecessary if this file is a symbolic link
from shared import *

# directories - main. change global_root with user_root as required
run_dir = os.path.join(global_root,'RunFolder')
user = getpass.getuser()
user_root = os.path.join(run_dir,user) #global_root
srf_default_dir = os.path.join(global_root, 'RupModel')
vel_mod_dir = os.path.join(global_root, 'VelocityModels')
vel_mod_subdirs = ['Cant','SI']
recipe_dir = os.path.join(bin_process_dir,"recipes")
v_mod_1d_dir = os.path.join(global_root,'VelocityModel','Mod-1D')
gmsa_dir = os.path.join(global_root,'groundMotionStationAnalysis')
stat_dir = os.path.join(global_root,'StationInfo')


def q_accept_custom_rupmodel():
    show_horizontal_line()
    print "Do you wish to use custom rupture files?"
    show_horizontal_line()
    return show_yes_no_question()
   


def q_custom_rupmodel_path():
    verified = False
    while not verified:
        rupture_path = raw_input("Enter path to custom Rupture Model (the parent of Srf/Stoch SRF directory): ")
        if not os.path.exists(rupture_path):
            print "Specified path doesn't exist. Try again"
            continue
        verified = True
    return rupture_path
      
        
       
def q_select_rupmodel_dir(srf_dir):
    show_horizontal_line()
    print "Select Rupture Model - Step 1."
    show_horizontal_line()
    srf_options = os.listdir(srf_dir)
    srf_options.sort()
    srf_selected = show_multiple_choice(srf_options)
    print srf_selected
    srf_selected_dir = os.path.join(srf_dir,srf_selected,"Srf")
    try:
        srf_file_options = os.listdir(srf_selected_dir)
    except OSError:
        print "No such Srf directory : %s" %srf_selected_dir
        sys.exit()

    params_src_path=os.path.join(srf_dir,srf_selected,params_src)
    if not os.path.exists(params_src_path):
        print "Error: %s doesn't exist" %(params_src_path)
        sys.exit()

    return srf_selected,srf_selected_dir,srf_file_options,params_src_path



def q_select_rupmodel(srf_selected_dir,srf_file_options):
    show_horizontal_line()
    print "Select Rupture Model - Step 2."
    show_horizontal_line()
    srf_file_options.sort()
    srf_files_selected = show_multiple_choice(srf_file_options,singular=False) #always in list
    
    print srf_files_selected
    srf_stoch_pairs=[]
    for srf_file_selected in srf_files_selected:
        srf_file_path = os.path.join(srf_selected_dir,srf_file_selected)
        stoch_selected_dir = os.path.abspath(os.path.join(srf_selected_dir,os.path.pardir,'Stoch'))
        if not os.path.isdir(stoch_selected_dir):
            print "Error: Corresponding Stoch directory is not found:\n%s" %stoch_selected_dir
            sys.exit()
        stoch_file_path = os.path.join(stoch_selected_dir,srf_file_selected).replace(".srf",".stoch")
        if not os.path.isfile(stoch_file_path):
            print "Error: Corresponding Stock file is not found:\n%s" %stoch_file_path
            sys.exit()

        print "Corresponding Stock file is also found:\n%s" %stoch_file_path
        srf_stoch_pairs.append((srf_file_path,stoch_file_path))

    return srf_files_selected, srf_stoch_pairs


#def q3():
#    show_horizontal_line()
#    print "Select HH "
#    show_horizontal_line()
#    HH_options = ['0.100','0.200','0.400','0.500']
#    HH = show_multiple_choice(HH_options)
#    print HH
#    return HH


def q_select_vel_model(vel_mod_dir):
    show_horizontal_line()
    print "Select one of available VelocityModels (from %s)" %vel_mod_dir
    show_horizontal_line()
    v_mod_ver_options = [] 
    for subdir in vel_mod_subdirs:
        vmodels=os.listdir(os.path.join(vel_mod_dir,subdir))
        v_mod_ver_options.extend([os.path.join(subdir,x) for x in vmodels]) 
    v_mod_ver_options.sort()
    v_mod_ver = show_multiple_choice(v_mod_ver_options)
    print v_mod_ver
    vel_mod_dir = os.path.join(vel_mod_dir,v_mod_ver)
    print vel_mod_dir
    params_vel_path = os.path.join(vel_mod_dir,params_vel)
    if not os.path.exists(params_vel_path):
        print "Error: %s doesn't exist" %(params_vel_path)
        sys.exit()

    return v_mod_ver,vel_mod_dir, params_vel_path

def q_select_stat_file(remove_fd=False):
    show_horizontal_line()
    print "Select one of available Station list (from %s)" %stat_dir
    show_horizontal_line()
    stat_files = glob.glob(os.path.join(stat_dir,'*.ll'))
    stat_files = [os.path.basename(x) for x in stat_files]
    if remove_fd:
        stat_files = [x for x in stat_files if not x.startswith('fd_')]
    
    stat_files.sort()
    stat_file = show_multiple_choice(stat_files)
    print stat_file
    stat_file_path = os.path.join(stat_dir,stat_file)
    print stat_file_path
    return stat_file_path


def q_get_run_name(HH,srf_selected,v_mod_ver,emod3d_version):
    #automatic generation of the run name (LP here only, HF and BB come later after declaration of HF and BB parameters). 
    userString=datetime.date.today().strftime("%y%m%d")   #additional string to customize (today's date for starters)
    hString='-h'+HH
    srfString=srf_selected #use full name of RupModel directory #srf_selected.split("_")[0]
    vModelString='VM'+str(v_mod_ver)
    vString='_EMODv'+emod3d_version
    run_name=(srfString+'_'+vModelString+hString+vString+'_'+userString).replace('.','p').replace('/','_')     #replace the decimal points with p
    # LPSim-2010Sept4_bev01_VMv1p64-h0p100_EMODv3p0p4_19May2016

    yes = confirm_name(run_name)
    return yes, run_name


def q_select_recipe(recipe_dir):
    recipes = os.listdir(recipe_dir)
    show_horizontal_line()
    print "Choose one of available recipes (from %s)" % recipe_dir
    show_horizontal_line()
    recipe_selected = show_multiple_choice(recipes)
    print recipe_selected
    recipe_selected_dir = os.path.join(recipe_dir,recipe_selected)
    return recipe_selected_dir


def q_final_confirm(run_name,recipe_selected_dir):
    show_horizontal_line(c="*")
    print "To be created: \n%s" %run_name
    print "Recipe to be copied from \n%s" %recipe_selected_dir
    show_horizontal_line(c="*")

    print "Do you wish to proceed?"
    return show_yes_no_question()


     
def action(sim_dir,recipe_selected_dir,run_name,version, global_root, user_root, run_dir, vel_mod_dir,srf_dir,srf_stoch_pairs,params_src_path,params_vel_path,stat_file_path):

    lf_sim_root_dir, hf_dir, bb_dir, figures_dir  = os.path.join(sim_dir,"LF"), os.path.join(sim_dir,"HF"), os.path.join(sim_dir,"BB"), os.path.join(sim_dir,"Figures")

    dir_list = [sim_dir, lf_sim_root_dir, hf_dir, bb_dir, figures_dir]
    if not os.path.isdir(user_root):
	dir_list.insert(0,user_root)
 
    verify_user_dirs(dir_list)

    for filename in glob.glob(os.path.join(recipe_selected_dir, '*.*')):
        shutil.copy(filename, sim_dir)

    shutil.copy(os.path.join(gmsa_dir,"parametersStation.py"),sim_dir)
    shutil.copy(os.path.join(gmsa_dir,"runPostProcessStation.ll"),sim_dir)
#    execute_cmd('ln -s %s/submit_emod3d.py %s'%(bin_process_dir,sim_dir))
    shutil.copy(os.path.join(bin_process_dir,"version"),sim_dir)
    shutil.copy(os.path.join(bin_process_dir,"submit_emod3d.sh"),sim_dir)
    shutil.copy(os.path.join(bin_process_dir,"submit_post_emod3d.sh"),sim_dir)
    shutil.copy(os.path.join(bin_process_dir,"submit_bb.sh"),sim_dir)

    srf_files, stoch_files = zip(*srf_stoch_pairs)
    with open(os.path.join(sim_dir,"params_base.py"),"w") as f:
        f.write("run_name='%s'\n" %run_name)
        f.write("version='%s'\n" %version)
        f.write("bin_process_ver='%s'\n" %bin_process_ver)

        f.write("global_root='%s'\n" %global_root)
        f.write("user_root='%s'\n" %user_root)
        f.write("run_dir='%s'\n"%run_dir)
        f.write("sim_dir='%s'\n"%sim_dir)
        f.write("lf_sim_root_dir='%s'\n"%lf_sim_root_dir)
        f.write("hf_dir='%s'\n"%hf_dir)
        f.write("bb_dir='%s'\n"%bb_dir)
        f.write("figures_dir='%s'\n"%figures_dir)
        f.write("srf_dir='%s'\n"%srf_dir)
        f.write("srf_files=%s\n"%str(list(srf_files)))
        f.write("hf_slips=%s\n"%str(list(stoch_files)))
        f.write("vel_mod_dir='%s'\n"%vel_mod_dir)
        f.write("v_mod_1d_dir='%s'\n"%v_mod_1d_dir)
        f.write("params_vel='%s'\n"%params_vel_path) #assumes one params_vel for all velocity models tested (not yet supported)
        f.write("params_src='%s'\n"%params_src_path) #assumes one params_src for all rupture models tested (supported)
        f.write("stat_file='%s'\n" %stat_file_path)
    

    print dir_list[0]
    set_permission(dir_list[0]) #if user_root is first time created, recursively set permission from there. otherwise, set permission from sim_dir

    print "Producing statcords and FD_STATLIST"
    sys.path.append(sim_dir)
    import statlist2gp
    fd_statcords, fd_statlist = statlist2gp.main()
    with open(os.path.join(sim_dir,"params_base.py"),"a") as f:
        f.write("stat_coords='%s'\n"%fd_statcords)
        f.write("FD_STATLIST='%s'\n"%fd_statlist)


def show_instruction(sim_dir):
    show_horizontal_line()
    print "Instructions"
    show_horizontal_line()
    print "    1.   cd %s" %sim_dir
    print "    2.   Edit params.py and run_emod3d.ll.template"
    print "    3.   ./submit_emod3d.sh"
    print "    4.   ./submit_post_emod3d.sh"
    print "    5.   (Linux) plot_and_ani.sh"
    print "    6.   install_bb.sh"
    print "    7.   ./submit_bb.sh"
    print "    8.   Edit parametersStation.py"
    print "    9.   llsubmit runPostProcessStation.ll"




def main():

    show_horizontal_line(c="*")
    print " "*37+"EMOD3D Job Preparation Ver."+bin_process_ver
    show_horizontal_line(c="*")
    
    yes = q_accept_custom_rupmodel()
    if yes:
        srf_dir = q_custom_rupmodel_path()
    else:
        srf_dir = srf_default_dir
    
    srf_selected,srf_selected_dir,srf_file_options,params_src_path = q_select_rupmodel_dir(srf_dir)
    srf_files_selected, srf_stoch_pairs = q_select_rupmodel(srf_selected_dir, srf_file_options)
#    HH = q3() ## HH is taken directly from params_vel.py
    v_mod_ver,vel_mod_dir_full, params_vel_path = q_select_vel_model(vel_mod_dir)
    
    execfile(params_vel_path,globals()) #import params_vel, to retrieve HH
    
    stat_file_path = q_select_stat_file(remove_fd=True) 
    
    yes, run_name = q_get_run_name(HH,srf_selected,v_mod_ver,emod3d_version)
    run_name = add_name_suffix(run_name,yes)
    
    recipe_selected_dir= q_select_recipe(recipe_dir)
    final_yes = q_final_confirm(run_name,recipe_selected_dir)
    if not final_yes:
        print "Installation exited"
        sys.exit()

    sim_dir = os.path.join(user_root,run_name)
    action(sim_dir,recipe_selected_dir,run_name,emod3d_version, global_root, user_root, run_dir, vel_mod_dir_full, srf_dir,srf_stoch_pairs,params_src_path,params_vel_path,stat_file_path)

    #add bin_process to PATH if it is not already there
    if not bin_process_dir in os.environ['PATH']:
        f=open("/home/%s/.bashrc" %user,'a')
        f.write("export PATH=$PATH:%s\n" %os.path.join(bin_process_path,bin_process_ver)) 

        f.close()
        print "PATH was updated"
    

    print "Installation completed"
    show_instruction(sim_dir)


if __name__ == '__main__':
    main()



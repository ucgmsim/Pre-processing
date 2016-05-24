import os
import sys
import datetime
import glob
import shutil

emod3d_version = '3.0.4'
bin_process_ver = 'devel'

# things that everyone doesn't have, eg. binaries are within here

#run_dir = os.path.dirname(os.path.realpath(__file__))
run_dir = os.path.abspath(os.path.curdir)
global_root = os.path.realpath(os.path.join(run_dir,os.path.pardir))
user_root = global_root
bin_process_dir = os.path.join(global_root, 'Pre-processing/bin_process',bin_process_ver)
# if bin_process_ver was not hardcoded and given as "stable" the following can workout the real version.
#bin_process_ver = os.readlink(bin_process_dir)
#bin_process_dir = os.path.realpath(os.path.join(bin_process_dir,os.path.pardir,bin_process_ver))

#print bin_process_ver
#print bin_process_dir


# directories - main. change global_root with user_root as required

srf_dir = os.path.join(user_root, 'RupModel')
vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel')
recipe_dir = os.path.join(bin_process_dir,"recipes")

def make_dirs(dir_list, reset = False):
    for dir_path in dir_list:
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        elif reset:
            # empty directory
            shutil.rmtree(dir_path)
            os.makedirs(dir_path)   


def user_select(options):
    try:
        selected_number = input("Enter the number you wish to select (1-%d):" %len(options))
    except NameError:
        print "Check your input."
        user_select(options)
    else:
        try:
            selected_number = int(selected_number)
        except ValueError:
            print "Input should be a number."
            user_select(options)
        else:
            try:
                return selected_number
            except IndexError:
                print "Input should be a number in (1-%d)" %len(options)
                user_select(options)


def show_multiple_choice(options):
    for i,option in enumerate(options):
        print "%2d. %s" %(i+1 , option)
    selected_number = user_select(options)
    return options[selected_number-1]

def show_yes_no_question():
    options = ["Yes","No"]
    for i, option in enumerate(options):
        print "%2d. %s" %(i+1 , option)
    selected_number = user_select(options)
    return (selected_number == 1 ) #True if selected Yes


def show_horizontal_line(c='=',length=100):
    print c*length

def q1(srf_dir):
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

    return srf_selected,srf_selected_dir,srf_file_options

def q2(srf_selected_dir,srf_file_options):
    show_horizontal_line()
    print "Select Rupture Model - Step 2."
    show_horizontal_line()
    srf_file_options.sort()
    srf_file_selected = show_multiple_choice(srf_file_options)
    print srf_file_selected
    srf_file_path = os.path.join(srf_selected_dir,srf_file_selected)
    return srf_file_selected,srf_file_path


def q3():
    show_horizontal_line()
    print "Select HH "
    show_horizontal_line()
    hh_options = ['0.100','0.500']
    hh = show_multiple_choice(hh_options)
    print hh
    return hh


def q4(vel_mod_dir):
    show_horizontal_line()
    print "Select one of available CanterburyVelocityModels (from %s)" %vel_mod_dir
    show_horizontal_line()
    
    v_mod_ver_options = os.listdir(vel_mod_dir)
    v_mod_ver_options.sort()
    v_mod_ver = show_multiple_choice(v_mod_ver_options)
    print v_mod_ver
    vel_mod_dir = os.path.join(vel_mod_dir,v_mod_ver)
    print vel_mod_dir
    return v_mod_ver,vel_mod_dir


def q5(hh,srf_selected,srf_file_selected,v_mod_ver,emod3d_version):
    #automatic generation of the run name (LP here only, HF and BB come later after declaration of HF and BB parameters). 
    userString=datetime.date.today().strftime("%d%B%Y")   #additional string to customize (today's date for starters)
    hString='-h'+hh
    srfString=srf_selected.split("_")[0]+"_"+os.path.splitext(srf_file_selected)[0]
    vModelString='VM'+str(v_mod_ver)
    vString='_EMODv'+emod3d_version
    run_name=(srfString+'_'+vModelString+hString+vString+'_'+userString).replace('.','p')     #replace the decimal points with p
    #run_name=('LPSim-'+srfString+'_'+vModelString+hString+vString+'_'+userString).replace('.','p')     #replace the decimal points with p
    # LPSim-2010Sept4_bev01_VMv1p64-h0p100_EMODv3p0p4_19May2016

    yes = happy_name(run_name)
    return yes, run_name

def happy_name(run_name):
    show_horizontal_line()
    print "Automated Run Name: ",
    print run_name
    show_horizontal_line()
    print "Do you wish to proceed?"
    return show_yes_no_question()

def q6(run_name,yes):
    new_run_name = run_name
    print "Yes? ",yes
    while not yes:
        userString = raw_input("Add more text (will be appended to the name above) ")
        userString=userString.replace(" ","_")
        new_run_name= run_name+"_"+userString
        yes = happy_name(new_run_name)
    return new_run_name

def q7(recipe_dir):
    recipes = os.listdir(recipe_dir)
    show_horizontal_line()
    print "Choose one of available recipes (from %s)" % recipe_dir
    show_horizontal_line()
    recipe_selected = show_multiple_choice(recipes)
    print recipe_selected
    recipe_selected_dir = os.path.join(recipe_dir,recipe_selected)
    return recipe_selected_dir


def q8(run_name,recipe_selected_dir):
    show_horizontal_line(c="*")
    print "To be created: \n%s" %run_name
    print "Recipe to be copied from \n%s" %recipe_selected_dir
    show_horizontal_line(c="*")

    print "Do you wish to proceed?"
    return show_yes_no_question()


def action(sim_dir,recipe_selected_dir,run_name,version, global_root,run_dir, vel_mod_dir,srf_dir,srf_file):
    make_dirs([sim_dir, os.path.join(sim_dir,"LF"), os.path.join(sim_dir,"HF"), os.path.join(sim_dir,"BB"), os.path.join(sim_dir,"Figures")])

    for filename in glob.glob(os.path.join(recipe_selected_dir, '*.*')):
        shutil.copy(filename, sim_dir)

    f=open(os.path.join(sim_dir,"params_base.py"),"w");
    f.write("run_name='%s'\n" %run_name)
    f.write("version='%s'\n" %version)
    f.write("global_root='%s'\n" %global_root)
    f.write("run_dir='%s'\n"%run_dir)
    f.write("sim_dir='%s'\n"%sim_dir)
    f.write("srf_dir='%s'\n"%srf_dir)
    f.write("srf_file='%s'\n"%srf_file)
    f.write("vel_mod_dir='%s'\n"%vel_mod_dir)
    f.close()

def show_instruction(sim_dir):
    show_horizontal_line()
    print "Instructions"
    show_horizontal_line()
    print "    1.   cd %s" %sim_dir
    print "    2.   Edit params.py"
    print "    3.   llsubmit run_emod3d.ll"
    print "    4.   llsubmit post_emod3d.ll"
    print "    5.   (Linux) bash plot_ts.sh"


def main():

    show_horizontal_line(c="*")
    print " "*40+"EMOD3D Job Preparation"
    show_horizontal_line(c="*")

    srf_selected,srf_selected_dir,srf_file_options = q1(srf_dir)
    srf_file_selected, srf_file = q2(srf_selected_dir, srf_file_options)
    hh = q3()
    v_mod_ver,vel_mod_dir_full = q4(vel_mod_dir)
    yes, run_name = q5(hh,srf_selected,srf_file_selected,v_mod_ver,emod3d_version)
    run_name = q6(run_name,yes)
    
    recipe_selected_dir= q7(recipe_dir)
    final_yes = q8(run_name,recipe_selected_dir)
    if not final_yes:
        print "Installation exited"
        sys.exit()

    sim_dir = os.path.join(run_dir,run_name)
    action(sim_dir,recipe_selected_dir,run_name,emod3d_version, global_root, run_dir, vel_mod_dir_full, srf_dir,srf_file)
    print "Installation completed"
    show_instruction(sim_dir)


if __name__ == '__main__':
    main()



import os
import sys
import datetime

emod3d_version = '3.0.4'

# things that everyone doesn't have, eg. binaries are within here

#run_dir = os.path.dirname(os.path.realpath(__file__))
run_dir = os.path.abspath(os.path.curdir)
global_root = os.path.realpath(os.path.join(run_dir,os.path.pardir))
user_root = global_root

bin_process_dir = os.path.join(global_root, 'Pre-processing/bin_process/stable')
print bin_process_dir
bin_process_ver = os.readlink(bin_process_dir)
bin_process_dir = os.path.realpath(os.path.join(bin_process_dir,os.path.pardir,bin_process_ver))

print bin_process_ver
print bin_process_dir


# directories - main. change global_root with user_root as required

srf_dir = os.path.join(user_root, 'RupModel')

v_mod_ver = 'v1.64'
vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel', v_mod_ver)
hh = '0.100'

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

def main():

#    if(len(sys.argv) <6):
#        print "Usage: %s srf_file hh velocity_model_path v_mod_ver userString....." %sys.argv[0]
#        print "   where"
#        print "     srf_file path:  eg. bev01.srf"
#        print "     velocity_model_path: dirname of rho3dfile.d, vp3dfile.p, vs3dfile.s (eg:...CanterburyVelocityModel/v1.64)
#        print "     hh: eg. 0.100"

#        sys.exit()
	
#    hh = sys.argv[1]
#    srf_file = sys.argv[2]
#    v_mod_ver = sys.argv[3]
#    version = sys.argv[4]
#    userString = '_'.join(sys.argv[5:])
    print "===================================="
    print "Select Rupture Model - Step 1."
    print "===================================="
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

    print "===================================="
    print "Select Rupture Model - Step 2."
    print "===================================="
    srf_file_options.sort()
    srf_file_selected = show_multiple_choice(srf_file_options)
    print srf_file_selected
    

    #automatic generation of the run name (LP here only, HF and BB come later after declaration of HF and BB parameters). 
    userString=datetime.date.today().strftime("%d%B%Y")   #additional string to customize (today's date for starters)
    hString='-h'+hh
    srfString=srf_selected.split("_")[0]+"_"+os.path.splitext(srf_file_selected)[0]
    vModelString='VM'+str(v_mod_ver)
    vString='_EMODv'+emod3d_version
    run_name=('LPSim-'+srfString+'_'+vModelString+hString+vString+'_'+userString).replace('.','p')     #replace the decimal points with p


    # LPSim-2010Sept4_bev01_VMv1p64-h0p100_EMODv3p0p4_19May2016
    print "===================================="
    print "Automated Run Name: ",
    print run_name
    print "===================================="
    print "Do you wish to proceed?"
    yes= show_yes_no_question()
    if not yes:
        userString = raw_input("Add more text (will be appended to the name above) ")
        run_name= run_name+"_"+userString
            
        print "===================================="
        print "Revised Run Name: ",
        print run_name
        print "===================================="
        

    

    recipe_dir = os.path.join(bin_process_dir,"recipes")
    print recipe_dir
    recipes = os.listdir(recipe_dir)

    print "===================================="
    print "Selected Recipes"
    print "===================================="
    recipe_selected = show_multiple_choice(recipes)
    print recipe_selected

    print "====================================="
    print "Directory %s created" %run_name
    print "Recipes from %s copied" %os.path.join(recipe_dir,recipe_selected)
    print "====================================="




if __name__ == '__main__':
    main()



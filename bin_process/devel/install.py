import os
import sys
import datetime

# things that everyone doesn't have, eg. binaries are within here
global_root = '/nesi/projects/nesi00213'


user_root = global_root

# works on Windows and POSIX paths
user_scratch = os.path.join(user_root, 'scratch', os.getenv('USER'))

# directories - main. change global_root with user_root as required
run_dir = os.path.join(user_root, 'RunFolder')
srf_dir = os.path.join(user_root, 'RupModel')
stat_dir = os.path.join(user_root, 'StationInfo')
vel_mod_params_dir = os.path.join(user_root, 'VelocityModel/ModelParams')
v_mod_ver = 'v1.64'
vel_mod_dir = os.path.join(global_root, 'CanterburyVelocityModel', v_mod_ver)
wcc_prog_dir = os.path.join(global_root, 'EMOD3D/WccFormat/bin/powerpc-AIX-nesi2-xlc')

# files
srf_file = os.path.join(srf_dir, '2010Sept4_m7pt1/Srf/bev01.srf')
print srf_file

stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_nz01-h0.100.statcords')    #TODO automate
hh = '0.100'
version = '3.0.4'

def main():

#    if(len(sys.argv) <6):
#        print "Usage: %s srf_file hh v_mod_ver version userString....." %sys.argv[0]
#        print "     where"
#        print "         srf_file path:  eg. bev01.srf"
#        print "         hh: eg. 0.100"
#        print "         v_mod_ver: v1.64
#        sys.exit()
	
#    hh = sys.argv[1]
#    srf_file = sys.argv[2]
#    v_mod_ver = sys.argv[3]
#    version = sys.argv[4]
#    userString = '_'.join(sys.argv[5:])


    script_dir = os.path.dirname(os.path.realpath(__file__))
    print script_dir
    recipe_dir = os.path.join(script_dir,"recipes")
    print recipe_dir
    recipes = os.listdir(recipe_dir)


    #automatic generation of the run name (LP here only, HF and BB come later after declaration of HF and BB parameters). 
    userString=datetime.date.today().strftime("%d%B%Y")   #additional string to customize (today's date for starters)
    hString='-h'+hh
    srfString=srf_file.split('RupModel/')[-1].split('_')[0]+'_'+srf_file.split('/')[-1].split('.')[0] 
    vModelString='VM'+str(v_mod_ver)
    vString='_EMODv'+version
    run_name=('LPSim-'+srfString+'_'+vModelString+hString+vString+'_'+userString).replace('.','p')     #replace the decimal points with p

    print run_name


if __name__ == '__main__':
    main()



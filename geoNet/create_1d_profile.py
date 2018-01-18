#!/usr/bin/python
"""create 1D profile"""
import os
import sys
import subprocess as sp
import argparse
import shutil
import glob

WORK_DIR = ".Multiple_Profiles"
OUTPUT_DIR = "Output"
OLD_ONED_DIR = "Output/Profiles"
ONED_DIR = "1D"
VM_DATA = "Data"
COORDINATE_TEXTFILE = "ProfileSites1D.txt"
DEPTH_PT_TEXTFILE = "ProfileDepthPoints1D.txt"
PROFILE = "GENERATE_MULTIPLE_PROFILES.txt"
DEPTH_PTS = [34, 0.025, 0.075, 0.125, 0.175, 0.275, 0.375, 0.575, 0.775, 0.975, 1.175, 1.375, 1.575, 1.775, 1.975,
             2.175, 2.375, 2.575, 2.775, 2.975, 3.175, 3.375, 3.575, 3.775, 3.975, 4.175, 4.375, 4.575, 4.775, 4.975,
             7.975, 11.975, 26.975, 38.975, 100.000]

nzvm_bin='NZVM'
VM_DATA_FULL_PATH=''

def mk_dir(curdir,mk_dir_name):
    """make a sub dir inside the current dir and cd to the sub dir"""
    work_path = os.path.join(curdir, mk_dir_name)
    #    work_path=os.path.join(loc,WORK_DIR)
    print"make", mk_dir_name
    if os.path.isdir(work_path):
        shutil.rmtree(work_path)
    os.makedirs(work_path)
    os.chdir(work_path)


def write_profile():
    """set up the structure of the profile"""
    with open(PROFILE, 'w') as f:
        f.write("CALL_TYPE=GENERATE_MULTIPLE_PROFILES\n")
        f.write("MODEL_VERSION=1.66\n")
        f.write("OUTPUT_DIR=%s\n" % OUTPUT_DIR)
        f.write("OUTPUT_TYPE=1D_SITE_RESPONSE\n")
        f.write("PROFILE_MIN_VS=0.00\n")
        f.write("TOPO_TYPE=SQUASHED\n")
        f.write("COORDINATES_TEXTFILE=%s\n" % COORDINATE_TEXTFILE)
        f.write("SPACING_TYPE=VARIABLE\n")
        f.write("PROFILE_DEPTHS_TEXTFILE=ProfileDepthPoints1D.txt\n")


def make_symbolic_link():
    """"make a symbolic link to Data"""
    if os.path.exists(VM_DATA):
        if os.path.islink(VM_DATA):
            pass
        else:
            print "Error: A directory %s already exists" % VM_DATA
            sys.exit()
    else:
        try:
            prog = sp.Popen(["find_config.sh"], stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
        except:
            print "Using user specified VM_DATA_FULL_PATH: %s" %VM_DATA_FULL_PATH
            os.symlink(VM_DATA_FULL_PATH,VM_DATA)
        else:
            out, err = prog.communicate()
            run_dir = os.path.dirname(out)
            os.symlink(os.path.join(run_dir, "VM/Velocity-Model/%s" % VM_DATA), VM_DATA)


def delete_dir(dir_name):
    """delete Multiple_Profiles/Output dir"""
    if os.path.exists(dir_name):
        print "Warning: %s exists. Deleting it" % dir_name
        shutil.rmtree(dir_name)


def generate_mutilprofile_txt(ll_file):
    """generate MultipleProfileParameters.txt"""
    with open(ll_file, 'r') as f:
        lines = f.readlines()
        num_lines = len(lines)
        with open(COORDINATE_TEXTFILE, 'w') as g:
            g.write("%d\n" % num_lines)
            for line in lines:
                s = ' '.join(filter(None, line.split(" ")))  # remove unnecessary spaces
                g.write("%s" % s)


def write_depth_file():
    """Write our depth point file"""
    with open(DEPTH_PT_TEXTFILE, 'w') as g:
        g.write("%d\n" % DEPTH_PTS[0])
        for line in DEPTH_PTS[1:]:
            g.write("%.3f\n" % line)


def exe_NZVM():
    """execute NZVM GENERATE_MULTIPLE_PROFILES.txt"""
    prog = sp.Popen([nzvm_bin, PROFILE], stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    out, err = prog.communicate()
    print out
    print err


def mv_oned_file(curdir):
    print "move one d"
    print "11111",os.getcwd()
    path = os.path.join(curdir,WORK_DIR,OLD_ONED_DIR)
    print(path)
    files = glob.glob1(path,'*.1d')
    for f in files:
        filepath = os.path.join(path,f)
        shutil.move(filepath, '.')


def event_generate_multiple_profiles(ll_file):
    """
        generates 1D profile
        ll_file: Absolute path to the .ll source file for creating 1D profile
        eg. home/seb56/20171218_Faketown_m6p3_201712181533/Stat/20171218_Faketown_m6p3/20171218_Faketown_m6p3.ll
    """
    curdir = os.path.abspath(os.curdir)  # get current directory

    mk_dir(curdir,WORK_DIR)

    write_profile()

    make_symbolic_link()

    delete_dir(OUTPUT_DIR)

    generate_mutilprofile_txt(ll_file)

    write_depth_file()

    exe_NZVM()

    mk_dir(curdir,ONED_DIR)

    mv_oned_file(curdir)

    os.chdir(curdir)

    delete_dir(WORK_DIR)



def check_args(fpath):
    """check if fname exists"""
    if os.path.isdir(fpath):
        print "Directory {} is present".format(fpath)
    elif os.path.isfile(fpath):
        print "File {} is present".format(fpath)
    else:
        sys.exit("File or directory {} does not exist".format(fpath))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ll_file', type=str,
                        help="Absolute path to the .ll source file for creating 1D profile\n eg. /home/seb56/20171218_Faketown_m6p3_201712181533/Stat/20171218_Faketown_m6p3/20171218_Faketown_m6p3.ll")
    parser.add_argument('NZVM',type=str, help="Full path of the NZVM executable\n eg. /home/seb56/code/vm/NZVM")
    parser.add_argument('VM_data',type=str, help="Full path to Velocity-Model Data directory. If missing, download https://quakecoresoft.canterbury.ac.nz/seisfinder/private/gmsim/vm_data_latest.tar.gz and extract or git clone https://github.com/ucgmsim/Velocity-Model.git")

    ll_file = parser.parse_args().ll_file
    nzvm_bin=parser.parse_args().NZVM
    VM_DATA_FULL_PATH=parser.parse_args().VM_data

    ll_file=os.path.abspath(ll_file)
    nzvm_bin=os.path.abspath(nzvm_bin)
    VM_DATA_FULL_PATH=os.path.abspath(VM_DATA_FULL_PATH)

    check_args(ll_file)
    check_args(nzvm_bin)
    check_args(VM_DATA_FULL_PATH)

    event_generate_multiple_profiles(ll_file)


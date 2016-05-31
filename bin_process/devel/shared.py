"""
Module which contains shared functions/values.

@date 8 April 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
"""

import os
import shutil

# reads a parameter from the parameters file (e3d.par)
# should not be necessary as you can just 'from params import *' (params.py)
def par_value (variable):
    result = ''
    par_handle = open('e3d.par', 'r')
    for line in par_handle:
        if line.startswith(variable + '='):
            # keep going and get the last result
            result = line
    par_handle.close()
    return ''.join(result.split('=')[1:]).rstrip('\n')

# returns a list of stations
# sample line in source file:
#   171.74765   -43.90236 ADCS
def get_stations(source_file, locations = False):
    stations = []
    station_lats = []
    station_lons = []
    with open(source_file, 'r') as sp:
        for line in sp.readlines():
            if line[0] not in  ['#', '%']:
                info = line.split()
                stations.append(info[2])
                if locations:
                    station_lons.append(info[0])
                    station_lats.append(info[1])
    if not locations:
        return stations
    return (stations, station_lats, station_lons)

# returns a dictionary of vrefs or vsites
# sample line in source file:
#SITE   VALUE
def get_vs(source_file):
    vs = {}
    with open(source_file, 'r') as sp:
        for line in sp.readlines():
            if line[0] not in ['#', '%']:
                info = line.split()
                vs[info[0]] = info[1]
    return vs


################# Verify Section ###################
# verify functions make sure script resources exist before continuing to run.
# it also creates output directories if not existing.
# very important to prevent (for example) empty variables. 'rm -r emptyvar/*' == 'rm -r /*'
# these functions prevent a script malfunctioning and causing code to run dangerously

# exception to throw, prevents scripts from continuing
class ResourceError(Exception):
    pass

# makes sure script file resources exist
def verify_files(file_list):
    for file_path in file_list:
        if not os.path.isfile(file_path):
            raise ResourceError('File not found: %s. Check params.py.' % (file_path))

# makes sure logfiles can be created, removes old ones
def verify_logfiles(logfile_list):
    for logfile in logfile_list:
        # reformat if just filename without path
        if os.path.dirname(logfile) == '':
            logfile = os.path.join(os.getcwd(), logfile)
        # is directory writable?
        if not os.access(os.path.dirname(logfile), os.W_OK):
            raise ResourceError('Can\'t write logfile: %s. Check directory permissions.'\
                    % (logfile))
        if os.path.exists(logfile):
            os.remove(logfile)

# makes sure required string are not empty
def verify_strings(string_list):
    for variable in string_list:
        if variable == '':
            raise ResourceError('Variable is empty: %s. Check params.py.' % (variable))

# makes sure list inputs contain values
def verify_lists(list_list):
    for req_list in list_list:
        if len(req_list) < 1:
            raise ResourceError('List doesn\'t contain any values: %s. Check params.py.' % (req_list))

# makes sure dirs which should already exist, do exist
def verify_dirs(dir_list):
    for dir_path in dir_list:
        if not os.path.isdir(dir_path):
            raise ResourceError('Directory doesn\'t exist: %s. Check params.py' % (dir_path))

# makes sure user dirs (ones that may be created if not existing) are ready
def verify_user_dirs(dir_list, reset = False):
    for dir_path in dir_list:
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)
        elif reset:
            # empty directory
            shutil.rmtree(dir_path)
            os.makedirs(dir_path)

# makes sure binary paths are valid binaries
def verify_binaries(bin_list):
    for bin_path in bin_list:
        if not os.path.isfile(bin_path):
            raise ResourceError('Binary not found: %s. Check params.py.' % (bin_path))
        if not os.access(bin_path, os.X_OK):
            raise ResourceError('Binary not executable: %s. Check file permissions.' % (bin_path))



 
#recursively sets permission. mode should be given in 0o777 format. eg. 0o750
def set_permission(dir_path,mode):

    print "Permission %s : %o" %(dir_path,mode)
    os.chmod(dir_path,mode)
    for root, dirs, files in os.walk(dir_path):
        for d in dirs:
            print "Permission %s : %o" %(os.path.join(root,d),mode)
            os.chmod(os.path.join(root,d),mode) 
        for f in files:
            print "Permission %s : %o" %(os.path.join(root,f),mode)
            os.chmod(os.path.join(root,f),mode)


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
                v = options[selected_number-1]
            except IndexError:
                print "Input should be a number in (1-%d)" %len(options)
                selected_number = user_select(options)
            return selected_number
         


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

def confirm_name(name):
    show_horizontal_line()
    print "Automated Name: ",
    print name
    show_horizontal_line()
    print "Do you wish to proceed?"
    return show_yes_no_question()

def add_name_suffix(name,yes):
    new_name = name
    print "Yes? ",yes
    while not yes:
        userString = raw_input("Add more text (will be appended to the name above) ")
        userString=userString.replace(" ","_")
        new_name= name+"_"+userString
        yes = confirm_name(new_name)
    return new_name



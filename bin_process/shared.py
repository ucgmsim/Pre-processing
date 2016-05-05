"""
Module which contains shared functions/values.
"""

import os

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

# makes sure user dirs (ones that may be created if not existing) are ready
def verify_user_dirs(dir_list):
    for dir_path in dir_list:
        if not os.path.isdir(dir_path):
            os.makedirs(dir_path)

# makes sure binary paths are valid binaries
def verify_binaries(bin_list):
    for bin_path in bin_list:
        if not os.path.isfile(bin_path):
            raise ResourceError('Binary not found: %s. Check params.py.' % (bin_path))
        if not os.access(bin_path, os.X_OK):
            raise ResourceError('Binary not executable: %s. Check file permissions.' % (bin_path))


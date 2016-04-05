#!/usr/bin/env python

import os
import sys
from subprocess import call
from time import strftime
from __future__ import print_function

# python time format string, imitate `date` command
time_format = '%a %b %d %H:%M:%S %Z %Y'

#set directory location of the  'merge_tsP3' exe
BPATH = os.path.expanduser('~rwg43') + '/Bin'

#
# FILEROOT is the "name" of the 3D FD run
# SIMDIR is the directory name where the 3D FD run output exists
# TSFILE is the file created here with all of the time slices merged together
#
# NPROC should be the total number of CPUs used for the 3D FD run
#

FILEROOT = '2011Feb22_m6pt2bev01_Cantv1.64'
SIMDIR = 'OutBin'
TSFILE = SIMDIR + '/' + FILEROOT + '_xyts.e3d'

NPROC = 2048

# end of input variables

OUTLOG = 'merge.out'
FILELIST = 'tmp.filelist'

if os.path.isfile(FILELIST):
    os.remove(FILELIST)
if os.path.exists(FILELIST):
    raise IOError('FILELIST is not a file or cannot be removed!')

# file number to look for
n = 0
# number of files found
nf = 0


# open found file list file for appending
try:
    list_handle = open(FILELIST, 'a')
except IOError:
    print('FILELIST cannot be opened to append data! ' + __file__, file = sys.stderr)
    raise

# add files to file list
while n < NPROC:
    file = SIMDIR + '/' + FILEROOT + '_xyts-' + str(n).zfill(5) + '.e3d'
    if os.path.exists(file):
        list_handle.write(file)
        nf += 1
    n += 1

list_handle.close()


# open log file to append standard output from subprocess
try:
    log_handle = open(OUTLOG, 'a')
except IOError:
    print('OUTLOG cannot be opened to append data! ' + __file__, file = sys.stderr)
    raise

# execute binary, pipe output to log
call([BPATH + '/merge_tsP3', 'filelist=' + FILELIST, 'outfile=' + TSFILE, 'nfiles=' + nf], stdout = log_handle)

log_handle.close()


print(__file__ + ' completed on ' + strftime(time_format))


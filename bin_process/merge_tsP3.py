#!/usr/bin/env python

from __future__ import print_function
import os
import sys
from subprocess import call
from time import strftime

from shared import par_value

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

FILEROOT = par_value('name')
SIMDIR = par_value('main_dump_dir')
TSFILE = SIMDIR + '/' + FILEROOT + '_xyts.e3d'

NPROC = par_value('nproc')

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
        list_handle.write(file + '\n')
        nf += 1
    n += 1

list_handle.close()

# open log file to write standard output from subprocess
try:
    log_handle = open(OUTLOG, 'w')
except IOError:
    print('OUTLOG cannot be opened to write stderr! ' + __file__, file = sys.stderr)
    raise

# execute binary, pipe output to log
code = call([BPATH + '/merge_tsP3', 'filelist=' + FILELIST, 'outfile=' + TSFILE, 'nfiles=' + str(nf)], stderr = log_handle)

log_handle.close()

if code != 0:
    print('MERGE FAILED')
print(__file__ + ' completed on ' + strftime(time_format))


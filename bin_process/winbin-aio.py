#!/usr/bin/env python

import os
#from shutil import move
import shutil
from subprocess import call
from glob import glob

from params import *

if not os.path.exists(VELDIR):
    os.makedirs(VELDIR)
if not os.path.isdir(VELDIR):
    raise IOError('Output directory is not a directory!')


#sample line from fd statlist: (warning, last line has no '\n')
#  171.74765   -43.90236 ADCS
# create list of values to iterate over
stat_handle = open(FD_STATLIST, 'r')
LONS = []
LATS = []
STATS = []
for line in stat_handle:
    lon_lat_stat = filter(None, line.rstrip('\n').split(' '))
    LONS.append(lon_lat_stat[0])
    LATS.append(lon_lat_stat[1])
    STATS.append(lon_lat_stat[2])

COMPS = ['080', '170', 'ver']

# -------- END OF USER DEFINED INPUT -------------

AZ_COLUMN = 8

for stat in STATS:
    print(stat)


IX = ['0', '1', '2']
FLIP = ['1', '1', '-1']

list_handle = open(FILELIST, 'w')

filepattern= MAIN_OUTPDIR + '/'+EXTENDED_RUN_NAME+ '_seis*.e3d'
#print filepattern
list_of_files = '\n'.join([file_path for file_path in glob(filepattern)])
#print listoffiles
list_handle.write(list_of_files)

list_handle.close()

for s_index, stat in enumerate(STATS):
    # never used??
    cdist = -999
    vsite = -999

    print(LONS[s_index] + ' ' + LATS[s_index] + ' ' + stat)

    statfile = VELDIR + '/' + stat

    for c_index, comp in enumerate(COMPS):
        cmd = [WCC_PROGDIR + '/fdbin2wcc', 'all_in_one=1', 'filelist=' + FILELIST,
                'ix=' + IX[c_index], 'scale=' + SCALE, 'flip=' + FLIP[c_index], 'tst=' + TSTRT,
                'stat=' + stat, 'comp=' + comp, 'outfile=' + statfile + '.' + comp,
                'nseis_comps=9', 'swap_bytes=0', 'tst2zero=0', 'outbin=0'] 
        print ' '.join(cmd)
        call(cmd)

    cmd = [WCC_PROGDIR + '/wcc_rotate','filein1=' + statfile + '.' + COMPS[0],'inbin1=0', 'filein2=' +  statfile + '.' + COMPS[1],'inbin2=0','fileout1=' + statfile + '.000', 'outbin1=0','fileout2=' + statfile+'.090','outbin2=0', 'rot=0.0']
    print ' '.join(cmd)
    call(cmd)

    
    os.remove('%s.%s'%(statfile,COMPS[0]))
    os.remove('%s.%s'%(statfile,COMPS[1]))


    
    


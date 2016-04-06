#!/usr/bin/env python

import os
from shutil import move
from subprocess import call

# -------- START OF USER DEFINED INPUT -------------

# Define location of EXE files
EXE_PATH = os.path.expanduser('~rwg43') + '/Bin'

# Define location of FD output files to use
SEISDIR = 'OutBin'

# Define folder to write velocity seismograms to
VELDIR = 'Vel'
OUTDIR = VELDIR
RUN = '2011Feb22_m6pt2bev01_Cantv1.64'

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
if not os.path.isdir(OUTDIR):
    raise IOError('Output directory is not a directory!')

FILELIST = 'fdb.filelist'

# Define 'start' time for time-axis (to ensure causality in filering) and scale (typ 1.0)
TSTRT = '-1.0'
SCALE = '1.0'

# Define the directory with the station information, and components
STAT_DIR = os.path.expanduser('~bab70') + '/StationInfo'
FD_STATLIST = STAT_DIR + '/fd_nz01-h0.100.ll'
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
list_handle.write('\n'.join([os.path.basename(file_path) for file_path in
        glob(SEISDIR + '/' + RUN + '_seis*.e3d')]))
list_handle.close()

for s_index, stat in enumerate(STATS):
    # never used??
    cdist = -999
    vsite = -999

    print(LONS[s_index] + ' ' + LATS[s_index] + ' ' + stat)

    for c_index, comp in enumerate(COMPS):
        call([EXE_PATH + '/fdbin2wcc', 'all_in_one=1', 'filelist=' + FILELIST,
                'ix=' + IX[c_index], 'scale=' + SCALE, 'flip=' + FLIP[c_index], 'tst=' + TSTRT,
                'stat=' + stat, 'comp=' + comp, 'outfile=' + stat + '.' + comp,
                'nseis_comps=9', 'swap_bytes=0', 'tst2zero=0', 'outbin=0'])

    call([EXE_PATH + '/wcc_rotate', 'filein1=' + stat + '.' + COMPS[0], 'inbin1=0',
            'filein2=' + stat + '.' + COMPS[[1], 'inbin2=0',
            'fileout1=' + stat + '.000', 'outbin1=0',
            'fileout2=' + stat + '.090', 'outbin2=0',
            'rot=0.0'])

    # shutil.move much simpler than os.rename
    move(stat + '.000', OUTDIR)
    move(stat + '.090', OUTDIR)
    move(stat + '.ver', OUTDIR)

    os.remove(stat + '.' + COMPS[0])
    os.remove(stat + '.' + COMPS[1])


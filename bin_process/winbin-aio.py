#!/usr/bin/env python

import os

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
TSTRT = -1.0
SCALE = 1.0

# Define the directory with the station information, and components
STAT_DIR = os.path.expanduser('~bab70') + '/StationInfo'
FD_STATLIST = STAT_DIR + '/fd_nz01-h0.100.ll'
#sample line from fd statlist: (take out last column only, warning, last line has no CR)
#  171.74765   -43.90236 ADCS
# create list of values to iterate over
STATS = 
#STATS = `gawk '{print $3;}' $FD_STATLIST `
COMPS = ( 080 170 ver )

# -------- END OF USER DEFINED INPUT -------------

set AZ_COLUMN = 8

echo $STATS


set IX = ( 0 1 2 )
set FLIP = ( 1 1 -1 )

ls ${SEISDIR}/${RUN}_seis*.e3d > $FILELIST

set s = 0
foreach stat ( $STATS )
@ s ++

set lon = `gawk -v stat=$stat '{if(substr($0,1,1)!="#" && $3==stat)x=$1;}END{print x;}' $FD_STATLIST `
set lat = `gawk -v stat=$stat '{if(substr($0,1,1)!="#" && $3==stat)x=$2;}END{print x;}' $FD_STATLIST `
set cdist = -999
set vsite = -999

echo $lon $lat $stat

set a = 0
foreach comp ( $COMPS )
@ a ++

${EXE_PATH}/fdbin2wcc all_in_one=1 filelist=$FILELIST \
	  ix=$IX[$a] scale=$SCALE flip=$FLIP[$a] tst=$TSTRT \
	  stat=$stat comp=$comp outfile=${stat}.${comp} \
	  nseis_comps=9 swap_bytes=0 tst2zero=0 outbin=0

end

${EXE_PATH}/wcc_rotate filein1=${stat}.$COMPS[1] inbin1=0 \
	   filein2=${stat}.$COMPS[2] inbin2=0 \
	   fileout1=${stat}.000 outbin1=0 \
           fileout2=${stat}.090 outbin2=0 \
	   rot=0.0

\mv ${stat}.000 ${stat}.090 ${stat}.ver $OUTDIR
\rm ${stat}.$COMPS[1] ${stat}.$COMPS[2]

end

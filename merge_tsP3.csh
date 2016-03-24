#!/bin/csh

#set directory location of the  'merge_tsP3' exe
set BPATH = ~rwg43/Bin

#
# FILEROOT is the "name" of the 3D FD run
# SIMDIR is the directory name where the 3D FD run output exists
# TSFILE is the file created here with all of the time slices merged together
#
# NPROC should be the total number of CPUs used for the 3D FD run
#

set FILEROOT = AlpineM7.9SthHypo_v1_SthIslandv1.01-h0.400_v3.04
set SIMDIR = OutBin
set TSFILE = ${SIMDIR}/${FILEROOT}_xyts.e3d

set NPROC = 2016

# end of input variables

set OUTLOG = merge.out
set FILELIST = tmp.filelist

\rm $FILELIST
set n = 0
set nf = 0
while( $n < $NPROC )

set file = `echo $n ${SIMDIR}/${FILEROOT} | gawk '{printf "%s_xyts-%.5d.e3d\n",$2,$1;}'`

if( -e ${file}) then

echo $file >> $FILELIST
@ nf ++

endif

@ n ++
end

$BPATH/merge_tsP3 filelist=$FILELIST outfile=$TSFILE nfiles=$nf >& $OUTLOG

echo "Done   " `date`

exit

#! /bin/csh

set BPATH = /hpc/home/rwg43/WccFormat/Progs

#
# FILEROOT is the "name" of the 3D FD run
# SIMDIR is the directory name where the 3D FD run output exists
# TSFILE is the file created by merge_ts
#
# DXTS is the decimation along the X-axis for the time slice output (same as dxts in e3d.par)
# DYTS is the decimation along the Y-axis for the time slice output
# DZTS is the decimation along the Z-axis for the time slice output (=1 for xy time slice)
# HH is the grid spacing used for the 3D FD run
#
# TS_INC is the increment for output time slices (=2 means every other slice is output)
# TS_START is the starting time slice (0 is first one)
# TS_TOTAL is the total number of slices to generate
#
# LONLAT_OUT =1 means output triplet is (lon,lat,amplitude) for point in the time slice
# GRIDFILE is the file containing the local (x,y,z) coordinates for this 3D run
#
# ABSMAX =1 vector magnitude of all 3 components is output into <fileroot>.0
#        =0 3 components are output respectively into <fileroot>.0, <fileroot>.1, <fileroot>.2
#

set FILEROOT = AlpineM7.9SthHypo_v1_SthIslandv1.01-h0.400_v3.04
set SIMDIR = ../OutBin
set TSFILE = $SIMDIR/${FILEROOT}_xyts.e3d

set HH = 0.400
set DXTS = 5
set DYTS = 5
set DZTS = 1

set ABSMAX = 1

set TS_INC = 1
set TS_START = 0
set TS_TOTAL = 3000

set LONLAT_OUT = 1
set GRIDFILE = /hpc/home/bab70/VelocityModel/SthIsland/ModelParams/gridout_sinz01-h${HH}

set OUTDIR = TSFiles
set OUTFILE = $OUTDIR/$FILEROOT

set SWAP_BYTES = 0
set SCALE = 1.0

\mkdir -p $OUTDIR

set tscnt = $TS_START
while($tscnt < $TS_TOTAL )
    @ tsnum = ( $tscnt * $TS_INC )

    set outf = `echo $OUTFILE $tscnt | gawk '{printf "%s_ts%.4d\n",$1,$2;}'`

    echo $outf

    $BPATH/ts2xyz infile=$TSFILE outfile=$outf swap_bytes=$SWAP_BYTES \
       gridfile=$GRIDFILE xyts=1 scale=$SCALE ts=$tsnum trv=0 \
       dxts=$DXTS dyts=$DYTS dzts=$DZTS absmax=$ABSMAX \
       read_header=1 outbin=1 lonlat=$LONLAT_OUT geoproj=1

    @ tscnt ++
end

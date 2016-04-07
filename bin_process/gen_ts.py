#! /usr/bin/env python2

import os
from subprocess import call

BPATH = os.path.expanduser('~rwg43') + '/WccFormat/Progs'

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

FILEROOT = '2011Feb22_m6pt2bev01_Cantv1.64'
SIMDIR = '../OutBin'
TSFILE = os.path.join(SIMDIR, FILEROOT + '_xyts.e3d')

HH = '0.100'
DXTS = '5'
DYTS = '5'
DZTS = '1'

ABSMAX = '1'

TS_INC = 1
TS_START = 0
TS_TOTAL = 400

LONLAT_OUT = '1'
GRIDFILE = os.path.expanduser('~bab70') + '/VelocityModel/ModelParams/gridout_nz01-h' + HH

OUTDIR = 'TSFiles'
OUTFILE = os.path.join(OUTDIR, FILEROOT)

SWAP_BYTES = '0'
SCALE = '1.0'

if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

for tscnt in range(TS_START, TS_TOTAL):
    tsnum = str(int(tscnt * TS_INC))

    outf = OUTFILE + '_ts' + str(tscnt).zfill(4)
    print(outf)

    call([os.path.join(BPATH, 'ts2xyz'), 'infile=' + TSFILE, 'outfile=' + outf, \
            'swap_bytes=' = SWAP_BYTES, 'gridfile=' + GRIDFILE, 'xyts=1', 'scale=' + SCALE, \
            'ts=' + tsnum, 'trv=0', 'dxts=' + DXTS, 'dyts=' + DYTS, 'dzts=' + DZTS, \
            'absmax=' + ABSMAX, 'read_header=1', 'outbin=1', 'lonlat=' + LONLAT_OUT, 'geoproj=1']


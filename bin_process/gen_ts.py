#! /usr/bin/env python2
"""
Converts OutBin/*_xyts.e3d files to TSlice/TSFiles/*_ts0000.0 which can then be plotted.

@author Viktor Polak, Sung Bae
@date 8 April 2016

Replaces gen_ts.csh. Converted to Python, uses globally set variables.

USAGE: execute from current directory. Input and output basenames are set in params.py.

ISSUES: 
"""

import os
from subprocess import call
import sys
import os.path
sys.path.append(os.path.abspath(os.path.curdir))
from params import *


if not os.path.exists(ts_out_dir):
    os.makedirs(ts_out_dir)



for tscnt in range(TS_START, TS_TOTAL):
    tsnum = str(int(tscnt * TS_INC))

    outf = ts_out_prefix + '_ts' + str(tscnt).zfill(4)
    print(outf)

    cmd = [os.path.join(wcc_prog_dir, 'ts2xyz'), 'infile=' + ts_file, 'outfile=' + outf, \
            'swap_bytes=' + swap_bytes, 'gridfile=' + GRIDFILE, 'xyts=1', 'scale=' + scale, \
            'ts=' + tsnum, 'trv=0', 'dxts=' + DXTS, 'dyts=' + DYTS, 'dzts=' + DZTS, \
            'absmax=' + ABSMAX, 'read_header=1', 'outbin=1', 'lonlat=' + lonlat_out, 'geoproj=1']
    print ' '.join(cmd)
    call(cmd)

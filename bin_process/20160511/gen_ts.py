#! /usr/bin/env python2
"""
Converts OutBin/*_xyts.e3d files to TSlice/TSFiles/*_ts0000.0 which can then be plotted.

@author Viktor Polak, Sung Bae
@date 8 April 2016

Replaces gen_ts.csh. Converted to Python, uses globally set variables.

USAGE: execute from current directory. Input and output basenames are set in params.py.

ISSUES: 
"""

from subprocess import call
import sys
import os.path
sys.path.append(os.path.abspath(os.path.curdir))
from shared import *
from params import *

verify_files([ts_file, GRIDFILE])
verify_strings([ts_start, ts_total, ts_inc, ts_out_prefix, swap_bytes, scale, \
        dx_ts, dy_ts, dz_ts, ABSMAX, lonlat_out])
verify_user_dirs([ts_out_dir])

print ts_start+" "+ts_total
ts_start = int(ts_start)
ts_total = int(ts_total)
ts_inc = int(ts_inc)

for tscnt in range(ts_start, ts_total):
    tsnum = str(int(tscnt * ts_inc))

    outf = ts_out_prefix + '_ts' + str(tscnt).zfill(4)
    print(outf)

    cmd = [os.path.join(wcc_prog_dir, 'ts2xyz'), 'infile=' + ts_file, 'outfile=' + outf, \
            'swap_bytes=' + swap_bytes, 'gridfile=' + GRIDFILE, 'xyts=1', 'scale=' + scale, \
            'ts=' + tsnum, 'trv=0', 'dxts=' + dx_ts, 'dyts=' + dy_ts, 'dzts=' + dz_ts, \
            'absmax=' + ABSMAX, 'read_header=1', 'outbin=1', 'lonlat=' + lonlat_out, 'geoproj=1']
    print ' '.join(cmd)
    call(cmd)

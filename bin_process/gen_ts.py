#! /usr/bin/env python2

import os
from subprocess import call

from params import *


if not os.path.exists(TS_OUTFILE_DIR):
    os.makedirs(TS_OUTFILE_DIR)



for tscnt in range(TS_START, TS_TOTAL):
    tsnum = str(int(tscnt * TS_INC))

    outf = TS_OUTFILE_PREFIX + '_ts' + str(tscnt).zfill(4)
    print(outf)

    cmd = [os.path.join(WCC_PROGDIR, 'ts2xyz'), 'infile=' + TSFILE, 'outfile=' + outf, \
            'swap_bytes=' + SWAP_BYTES, 'gridfile=' + GRIDFILE, 'xyts=1', 'scale=' + SCALE, \
            'ts=' + tsnum, 'trv=0', 'dxts=' + DXTS, 'dyts=' + DYTS, 'dzts=' + DZTS, \
            'absmax=' + ABSMAX, 'read_header=1', 'outbin=1', 'lonlat=' + LONLAT_OUT, 'geoproj=1']
    print ' '.join(cmd)
    call(cmd)

#!/usr/bin/env python2
"""
Generates LON,LAT,VEL IEEE754 32bit (single precision) binary file.
Combines multiple TimeSlice files into a single one.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 June 2016

USAGE: 

ISSUES: currently assumes ABSMAX TS input.
    TODO: modify as required.
"""

from glob import glob
from math import ceil, sqrt

import numpy as np

from params import nx, ny, dx_ts, dy_ts

xs = nx // dx_ts
ys = ny // dy_ts
points = xs * ys

# input
tsfiles = glob('TSFiles/*.0')

# np.memmap cannot open more than 1000 files.
# work on smaller sets to reduce memory usage.
sset = 1000
tsmax = np.zeros(shape = (points))
for fset in xrange(int(ceil(len(tsfiles) / float(sset)))):
    tsmax = np.maximum(tsmax, np.max(np.array([np.fromfile(f, dtype = '3>f')[:, -1] \
            for f in tsfiles[fset * sset:(fset + 1) * sset]], dtype = 'f'), axis = 0))

# replace latitude and longitude
tsmaxgrid = np.fromfile(tsfiles[0], dtype = '>f')
tsmaxgrid[2::3] = tsmax

# output
tsmaxgrid.astype(np.float32).tofile('statgrid_max.bin', format = '>f4')

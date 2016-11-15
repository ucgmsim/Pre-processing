#!/usr/bin/env python2
"""
Generates a time series from a series of time slices.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 November 2016

USAGE: 

ISSUES: currently assumes ABSMAX TS input.
    TODO: modify as required.
"""

from glob import glob
from math import ceil, sqrt

import numpy as np

from params import nx, ny, dx_ts, dy_ts, nt, dt, dt_ts
from shared_ts import vel2acc

dt = float(dt) * int(dt_ts)
xs = int(nx) // int(dx_ts)
ys = int(ny) // int(dy_ts)
points = xs * ys

# input
tsfiles = sorted(glob('TSFiles/*.0'))

# find points of interest
use_points = []
f = np.fromfile(tsfiles[0], dtype = '3>f')
for i in xrange(points):
    lon, lat = f[i, :2]
    if 170.18415 < lon < 170.20025 and -43.3846666 > lat > -43.3961139:
        print('found point %d at %s, %s' % (i, lat, lon))
        use_points.append(i)

for p in use_points:
    ts = np.array([np.fromfile(f, dtype = '3>f')[p, -1] \
            for f in tsfiles])
    ts = vel2acc(ts, dt)
    # output
    ts.astype(np.float32).tofile('ts%d_le.bin' % (p))
    np.savetxt('ts%d.txt' % (p), ts, fmt = '%f', delimiter = ' ')

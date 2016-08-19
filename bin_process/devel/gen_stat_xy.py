#!/usr/bin/env python2

import sys

import numpy as np

from params import *
from tools import ll2gp, InputError

# have to know which point the station is at
lat = -43.5928
lon = 172.3811

try:
    x, y = ll2gp(lat, lon, float(MODEL_LAT), float(MODEL_LON), \
            float(MODEL_ROT), int(nx), int(ny), float(hh), \
            int(dx_ts), int(dy_ts), decimated = True)
except InputError:
    print('Station coordinates outside simulation domain')
    raise
# position within timeslice file which this x, y appears
pos = y * int(nx) / int(dy_ts) + x

# seismo data in timeslices, these have to be present for plotting anyway
# reduces dependencies of files, can plot stations not in statfile
ts = [ts_out_prefix + '_ts%.4d.0' % (i) \
        for i in xrange(int(ts_total))]
ds = [np.memmap(filename = ts[i], dtype = '3>f', mode = 'r')[pos][-1] \
        for i in xrange(int(ts_total))][:150]
v_max = max(np.abs(ds))


###
### GENERATE XY
###
xfac = 0.2/len(ds)
yfac = 0.1/v_max
with open('station.xy', 'w') as sp:
    for i, value in enumerate(ds):
        lat2 = lat + value * yfac
        lon2 = lon + i * xfac
        sp.write('%f %f\n' % (lon2, lat2))

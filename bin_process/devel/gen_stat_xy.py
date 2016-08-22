#!/usr/bin/env python2

from math import sin, cos, radians
import sys

import numpy as np

from params import *
from tools import *

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
        for i in xrange(int(ts_total))]
v_max = max(np.abs(ds))


###
### GENERATE XY
###
xfac = 35./len(ds)
yfac = 15./v_max
lat0 = lat + 0.01
lon0 = lon + 0.01
xazim = 90
if xazim % 360 - 180 < 0:
    yazim = xazim - 90
else:
    yazim = xazim + 90
lls = ['%f %f\n' % (lon0, lat0)]
for i, value in enumerate(ds):
    dy = value * yfac
    dx = i * xfac
    lat1, lon1 = ll_shift(lat0, lon0, dx, xazim)
    lat1, lon1 = ll_shift(lat1, lon1, dy, yazim)
    lls.append('%f %f\n' % (lon1, lat1))

with open('gmt-seismo.xy', 'w') as sp:
    sp.write(''.join(lls))

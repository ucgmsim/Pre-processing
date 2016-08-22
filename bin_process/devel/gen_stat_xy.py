#!/usr/bin/env python2

from math import sin, cos
import os
import sys

import numpy as np

from params import *
from tools import *

gmt_seis = 'gmt-seismo.xy'

if os.path.exists(gmt_seis):
    os.remove(gmt_seis)

# seismo data in timeslices, these have to be present for plotting anyway
# reduces dependencies of files, can plot stations not in statfile
ts = [ts_out_prefix + '_ts%.4d.0' % (i) \
        for i in xrange(int(ts_total))]

###
### read values, calculate maximums
###
all_ds = []
max_ds = 0
for desc in plot_seismo:
    lat, lon, offset, oazim, xazim, tlen, yamp = desc

    try:
        x, y = ll2gp(lat, lon, float(MODEL_LAT), float(MODEL_LON), \
                float(MODEL_ROT), int(nx), int(ny), float(hh), \
                int(dx_ts), int(dy_ts), decimated = True)
    except InputError:
        print('Station coordinates outside simulation domain')
        raise
    # position within timeslice file which this x, y appears
    pos = y * int(nx) / int(dy_ts) + x

    ds = [np.memmap(filename = ts[i], dtype = '3>f', mode = 'r')[pos][-1] \
            for i in xrange(int(ts_total))]

    all_ds.append(ds)
    v_max = max(np.abs(ds))
    max_ds = max(max_ds, v_max)

# y scaling factor based on overall max (km)
yfac = float(yamp)/max_ds

###
### generate XY data based on read values and max / other params
###
for i, desc in enumerate(plot_seismo):
    lat, lon, offset, oazim, xazim, tlen, yamp = desc

    # scaling factor to produce wanted x, y lengths (km)
    xfac = float(tlen)/len(all_ds[i])

    # offset of beginning of seismogram
    lat0, lon0 = ll_shift(lat, lon, offset, oazim)

    # rotate as wanted (do not visually flip seismo)
    if xazim % 360 - 180 < 0:
        yazim = xazim - 90
    else:
        yazim = xazim + 90

    lls = ['%f %f\n' % (lon0, lat0)]
    for i, value in enumerate(all_ds[i]):
        dy = value * yfac
        dx = i * xfac
        # find next point using distance / bearing
        lat1, lon1 = ll_shift(lat0, lon0, dx, xazim)
        lat1, lon1 = ll_shift(lat1, lon1, dy, yazim)
        lls.append('%f %f\n' % (lon1, lat1))

    with open(gmt_seis, 'a') as sp:
        sp.write('> station at %f %f\n' % (lon, lat))
        sp.write(''.join(lls))


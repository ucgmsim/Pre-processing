#!/usr/bin/env python2

from math import sin, cos
import os
import sys

import numpy as np

from params import *
from shared_ts import *
from tools import *

# output file
gmt_seis = 'gmt-seismo.xy'

if os.path.exists(gmt_seis):
    os.remove(gmt_seis)

# common properties from params.py
xazim, tlen, yamp = plot_seismo_params

###
### read values, calculate maximums
###
all_ds = []
max_ds = 0
xy = open('xy.ll', 'r')
for line in xy:
    lat, lon, stat = line.split()
    lat = float(lat)
    lon = float(lon)
    data = 'velBB/%s.000' % (stat)

    ds = read_ascii(data)

    all_ds.append(ds)
    # PGVs
    v_max = max(np.abs(ds))
    max_ds = max(max_ds, v_max)
xy.close()

# y scaling factor based on overall max (km)
yfac = float(yamp) / max_ds
xfac = float(tlen) / max([len(ts) for ts in all_ds])

###
### generate XY data based on read values and max / other params
### simple version where seismo line is extended
###
xy = open('xy.ll', 'r')
for i, line in enumerate(xy):
    lon, lat, name = line.split()
    lon = float(lon)
    lat = float(lat)
    offset, oazim = 0, 0

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
xy.close()

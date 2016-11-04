#!/usr/bin/env python2

from os import path

from tools import *
from shared import *
from params import *

verify_strings([MODEL_LAT, MODEL_LON, MODEL_ROT, hh, nx, ny])
verify_dirs([stat_dir])

# velocity model parameters
nx = int(nx)
ny = int(ny)
mlat = float(MODEL_LAT)
mlon = float(MODEL_LON)
mrot = float(MODEL_ROT)
hh = float(hh)

# arbitrary longlat station input
ll_in = stat_file
# where to save gridpoint and longlat station files
gp_out = path.join(stat_dir, 'fd_nz01-h%s.statcords' % (hh))
ll_out = path.join(stat_dir, 'fd_nz01-h%s.ll' % (hh))

# shortcut functions for conversion
def get_gp(lon, lat):
    return ll2gp(lat, lon, mlat, mlon, mrot, nx, ny, hh)
def get_ll(x, y):
    return gp2ll(x, y, mlat, mlon, mrot, nx, ny, hh)

# retrieve in station names, latitudes and longitudes
sname, slat, slon = get_stations(ll_in, locations = True)
slon = map(float, slon)
slat = map(float, slat)

# store gridpoints and names if unique position
sxy = []
suname = []
for i in xrange(len(sname)):
    xy = get_gp(slon[i], slat[i])
    if xy not in sxy:
        sxy.append(xy)
        suname.append(sname[i])
    else:
        print('Duplicate Station Ignored: %s' % (sname[i]))

# create gp file
with open(gp_out, 'w') as gpf:
    # file starts with number of entries
    gpf.write('%d\n' % (len(sxy)))
    # x, y, z, name
    for i, xy in enumerate(sxy):
        gpf.write('%5d %5d %5d %s\n' % (xy[0], xy[1], 1, suname[i]))

# creat ll file
with open(ll_out, 'w') as llf:
    # lon, lat, name
    for i, xy in enumerate(sxy):
        lon, lat = get_ll(xy[0], xy[1])
        llf.write('%11.5f %11.5f %s\n' % (lon, lat, suname[i]))

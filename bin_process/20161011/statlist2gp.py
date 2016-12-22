#!/usr/bin/env python2

from os import path

from tools import *
from shared import *
import sys
import os
sys.path.append(os.path.abspath(os.path.curdir))

try:
    from params import *
except ImportError:
    print "params.py is not available"
    print "Try params_base.py"
    from params_base import *



# shortcut functions for conversion
def get_gp(lon, lat,mlat,mlon,mrot,nx,ny,hh):
    return ll2gp(lat, lon, mlat, mlon, mrot, nx, ny, hh)
def get_ll(x, y,mlat,mlon,mrot,nx,ny,hh):
    return gp2ll(x, y, mlat, mlon, mrot, nx, ny, hh)


def main():
    
    if params_vel:
        execfile(params_vel,globals())
        verify_strings([ORIGIN_LAT,ORIGIN_LON,ORIGIN_ROT, HH, NX, NY,CODE,sim_dir])
        MODEL_LAT=ORIGIN_LAT
        MODEL_LON=ORIGIN_LON
        MODEL_ROT=ORIGIN_ROT
        hh = HH
        nx = NX
        ny = NY
        code = CODE
    else:
        verify_strings([MODEL_LAT, MODEL_LON, MODEL_ROT, hh, nx, ny,sim_dir])

    outpath = sim_dir
    filename = 'fd_%s01-h%s'%(code,hh)

    # arbitrary longlat station input
    ll_in = stat_file
    # where to save gridpoint and longlat station files
    gp_out = path.join(outpath, '%s.statcords' % (filename))
    ll_out = path.join(outpath, '%s.ll' % (filename))


    print "From: %s" %stat_file
    print "To:"
    print "  %s" %gp_out
    print "  %s" %ll_out

    # velocity model parameters
    nx = int(nx)
    ny = int(ny)
    mlat = float(MODEL_LAT)
    mlon = float(MODEL_LON)
    mrot = float(MODEL_ROT)
    hh = float(hh)

    # retrieve in station names, latitudes and longitudes
    sname, slat, slon = get_stations(ll_in, locations = True)
    slon = map(float, slon)
    slat = map(float, slat)

    # store gridpoints and names if unique position
    sxy = []
    suname = []
    for i in xrange(len(sname)):
        try:
            xy = get_gp(slon[i], slat[i],mlat,mlon,mrot,nx,ny,hh)
        except InputError:
            continue
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
            lon, lat = get_ll(xy[0], xy[1],mlat,mlon,mrot,nx,ny,hh)
            llf.write('%11.5f %11.5f %s\n' % (lon, lat, suname[i]))
    
    return gp_out, ll_out

if __name__ == "__main__":
    main()


"""
Various tools which may be needed in various processes.
"""

from math import sin, asin, cos, atan2, degrees, radians
from subprocess import Popen, PIPE

R_EARTH = 6378.139
# ideally implemented in python
ll2xy_bin = '/nesi/projects/nesi00213/tools/ll2xy'

class InputError(Exception):
    pass

def ll2gp(lat, lon, mlat, mlon, rot, nx, ny, hh, \
        dx = 1, dy = 1, decimated = False):
    """
    Converts latitude/longitude to a gridpoint position.
    Three main modes of operation:
    1: No dx, dy (= 1): gridpoint
    2: dx or dy != 1: closest gridpoint considering dx/dy
    3: decimated: gridpoint number if only decimated points existed
    """
    # derived parameters
    xlen = nx * hh
    ylen = ny * hh
    # where the x plane points
    xazim = (rot + 90) % 360

    # run binary, get output
    # output is displacement (x, y) from center, in kilometres
    p_conv = Popen([ll2xy_bin, 'mlat=%.5f' % (mlat), 'mlon=%.5f' % (mlon), \
            'geoproj=1', 'center_origin=1', 'h=%.5f' % (hh), \
            'xazim=%.5f' % (xazim), 'xlen=%f' % (xlen), 'ylen=%f' % (ylen)], \
            stdin = PIPE, stdout = PIPE)
    stdout = p_conv.communicate('%.15f %.15f' % (lon, lat))[0]
    x, y = map(float, stdout.split())

    # convert displacement to grid points
    # first make the distance relative to top corner
    max_x = (nx ) * hh #nx-1 if concerned of consistent behaviour handling the maximum (eg. 1399 -> 1400 is not allowed)
    max_y = (ny ) * hh
    x += (max_x) * 0.5
    y += (max_y) * 0.5
    # then convert back to grid spacing
    x /= hh
    y /= hh
    # gridpoints are discrete
    x = int(round(x))
    y = int(round(y))

    print "Raw XY coords: (x,y)=%d,%d" %(x,y)
    # nx values range from 0 -> nx - 1
    if not (-1 < x < nx) or not (-1 < y < ny):
        raise InputError('Input outside simulation domain.')

    if decimated:
        x //= dx
        y //= dy
    else:
        # closest gridpoint considering decimation
        x -= x % dx
        y -= y % dy

    return x, y

def ll_shift(lat, lon, distance, bearing):
    """
    Shift lat/long by distance at bearing.
    """
    # formula is for radian values
    lat, lon, bearing = map(radians, [lat, lon, bearing])

    shift = distance / R_EARTH
    lat2 = asin(sin(lat) * cos(shift) \
            + cos(lat) * sin(shift) * cos(bearing))
    lon2 = lon + atan2(sin(bearing) * sin(shift) * cos(lat), \
            cos(shift) - sin(lat) * sin(lat2))

    return degrees(lat2), degrees(lon2)

"""
Various tools which may be needed in various processes.
"""

from subprocess import Popen, PIPE

# ideally implemented in python
ll2xy_bin = '/home/vap30/bin/ll2xy'

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
    max_x = (nx - 1) * hh
    max_y = (ny - 1) * hh
    x += max_x * 0.5
    y += max_y * 0.5
    # then convert back to grid spacing
    x /= hh
    y /= hh
    # gridpoints are discrete
    x = int(round(x))
    y = int(round(y))

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

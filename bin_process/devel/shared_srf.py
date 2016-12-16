"""
Common SRF functions.

SRF format:
https://scec.usc.edu/scecpedia/Standard_Rupture_Format
"""

from math import ceil, cos, radians
from subprocess import Popen, PIPE

import numpy as np

from tools import *

# binary paths
srf2xyz = '/nesi/projects/nesi00213/tools/srf2xyz'
# assumption that all srf files contain 6 values per line
VPL = 6.

def get_nseg(srf):
    """
    Returns number of segments in SRF file.
    srf: filepath to srf
    """
    with open(srf, 'r') as sf:
        sf.readline()
        nseg = int(sf.readline().split()[1])
    return nseg

def read_header(sf):
    """
    Parse header information.
    sf: open srf file at position 0
    """
    version = float(sf.readline())
    nseg = int(sf.readline().split()[1])
    planes = []
    for _ in xrange(nseg):
        # works for version 1.0 and 2.0
        elon, elat, nstk, ndip, ln, wid = sf.readline().split()
        stk, dip, dtop, shyp, dhyp = sf.readline().split()
        # store as correct format
        planes.append((float(elon), float(elat), int(nstk), int(ndip), \
                float(ln), float(wid), float(stk), float(dip), \
                float(dtop), float(shyp), float(dhyp)))
    return planes

def skip_points(sf, np):
    """
    Skips wanted number of points entries in SRF.
    sf: open srf file at the start of a point
    np: number of points to read past
    """
    for _ in xrange(np):
        # header 1 not important
        sf.readline()
        # header 2 contains number of values which determines lines
        values = sum(map(int, sf.readline().split()[2::2]))
        for _ in xrange(int(ceil(values / VPL))):
            sf.readline()

def get_lonlat(sf, value = None):
    """
    Returns only the longitude, latitude of a point.
    sf: open file at start of point
    value: also retrieve value
    end: sf at start of next point
    """
    # header 1 contains:
    # LON, LAT, DEP, STK, DIP, AREA, TINIT, DT, VS (v2.0), DEN (v2.0)
    h1 = sf.readline().split()
    # header 2 contains:
    # RAKE, SLIP1, NT1, SLIP2, NT2, SLIP3, NT3
    h2 = sf.readline().split()

    # always returning lon, lat
    lon, lat = map(float, h1[:2])

    if value == 'tinit':
        value = h1[6]

    # skip rest of point data
    values = sum(map(int, h2[2::2]))
    for _ in xrange(int(ceil(values / VPL))):
        sf.readline()

    if value == None:
        return lon, lat
    return lon, lat, value

def get_bounds(srf, seg = -1):
    """
    Return corners of segments.
    srf: srf source
    nseg: which segment (-1 for all)
    """
    bounds = []
    with open(srf, 'r') as sf:
        # metadata
        planes = read_header(sf)
        points = int(sf.readline().split()[1])

        # each plane has a separate set of corners
        for n, plane in enumerate(planes):
            plane_bounds = []
            nstk, ndip = plane[2:4]
            # set of points starts at corner
            plane_bounds.append(get_lonlat(sf))
            # travel along strike, read last value
            skip_points(sf, nstk - 2)
            plane_bounds.append(get_lonlat(sf))
            # go to start of strike at bottom of dip
            skip_points(sf, (ndip - 2) * nstk)
            plane_bounds.append(get_lonlat(sf))
            # travel along strike at bottom of dip
            skip_points(sf, nstk - 2)
            plane_bounds.insert(2, get_lonlat(sf))
            # store plane bounds or return if only 1 wanted
            if n == seg:
                return plane_bounds
            bounds.append(plane_bounds)
    return bounds

def get_hypo(srf):
    """
    Return hypocentre.
    srf: srf source
    """
    with open(srf, 'r') as sf:
        # metadata
        plane = read_header(sf)[0]
        lon, lat = plane[:2]
        strike, dip = plane[6:8]
        # plane 9 11?? dum dum dum....
        shyp, dhyp = plane[9:11]

        # move along strike for shyp km
        lat, lon = ll_shift(lat, lon, shyp, strike)
        # move along dip for dhyp km
        flat_dhyp = dhyp * cos(radians(dip))
        lat, lon = ll_shift(lat, lon, flat_dhyp, strike + 90)

        return lon, lat

def srf2llv(srf, seg = '-1', type = 'slip', depth = False):
    """
    Get longitude, latitude, depth (optional) and value of 'type'
    srf: filepath of SRF file
    seg: which segments to read (-1 for all)
    type: which parameter to read
    depth: whether to also include depth at point
    """
    proc = Popen([srf2xyz, 'calc_xy=0', 'lonlatdep=1', \
            'dump_slip=0', 'infile=%s' % (srf), \
            'type=%s' % (type), 'nseg=%d' % (seg)], stdout = PIPE)
    out, err = proc.communicate()
    code = proc.wait()
    # process output
    llv = np.fromstring(out, dtype = 'f4', sep = ' ')

    # create a slice filter if depth not wanted
    # longitude, latitude, depth, value
    if not depth:
        mask = np.array([True, True, False, True])
    else:
        mask = np.array([True, True, True, True])

    # output from srf2xyz is 4 columns wide
    return np.reshape(llv, (len(llv) // 4, 4))[:, mask]

def get_tinit(srf, seg = -1):
    """
    Return lon, lat, tinit for subfaults.
    Should be part of srf2llv in the future.
    srf: srf source
    nseg: which segment (-1 for all)
    """
    with open(srf, 'r') as sf:
        # metadata
        planes = read_header(sf)
        points = int(sf.readline().split()[1])

        # storage
        tinit = []

        # each plane has a separate set of subfaults
        for n, plane in enumerate(planes):
            nstk, ndip = plane[2:4]
            if seg >= 0 and seg != n:
                skip_points(sf, nstk * ndip)
                continue

            plane_tinit = np.zeros((nstk * ndip, 3))
            for i in xrange(nstk * ndip):
                plane_tinit[i] = get_lonlat(sf, value = 'tinit')
            tinit.append(plane_tinit)

            if n == seg:
                break
    return tinit

def srf_dxy(srf):
    """
    Retrieve SRF dx and dy.
    Assumes all planes have same dx, dy.
    srf: SRF file path to read from
    """
    with open(srf, 'r') as sf:
        # version
        sf.readline()
        # planes definition
        sf.readline()
        # first plane
        elon, elat, nstk, ndip, length, width = sf.readline().split()
    return float('%.2f' % (float(length) / int(nstk))), \
            float('%.2f' % (float(width) / int(ndip)))

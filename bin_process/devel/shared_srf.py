"""
Common SRF functions.

SRF format:
https://scec.usc.edu/scecpedia/Standard_Rupture_Format
"""

from math import ceil, cos, radians

from tools import *

# assumption that all srf files contain 6 values per line
VPL = 6.

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

def get_lonlat(sf):
    """
    Returns only the longitude, latitude of a point.
    sf: open file at start of point
    end: sf at start of next point
    """
    lon, lat = map(float, sf.readline().split()[:2])
    # skip rest of point data
    values = sum(map(int, sf.readline().split()[2::2]))
    for _ in xrange(int(ceil(values / VPL))):
        sf.readline()
    return lon, lat

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
        for plane in planes:
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
            # store plane bounds
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


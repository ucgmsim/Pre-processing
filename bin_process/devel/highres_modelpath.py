#!/usr/bin/env python2
"""
Reads a model_params file to get simulation domain corners.
Generates a path with many more points inbetween corners.
Allows gmt grdmasks to behave correctly (not curve incorrectly).

@author: Viktor Polak
@contact: viktor.polak@canterbury.ac.nz
@date: 16 November 2016

ISSUES: slow but very simple algorithm
        (not computationally intensive for sane min_edge_points)
"""

from tools import ll_mid

def path_from_corners(corners = None, output = 'sim.modelpath_hr', min_edge_points = 100):
    """
    corners: python list (4 by 2) containing (lon, lat) in order
        otherwise take from velocity model
    output: where to store path of (lon, lat) values
    min_edge_points: at least this many points wanted along edges
    """
    # input data using velocity model
    if corners == None:
        # don't fail importing if not needed
        from params import *
        # load model_params
        with open(vel_mod_params, 'r') as mp:
            lines = mp.readlines()
        # find corners using tags
        for tag in ['c1=', 'c2=', 'c3=', 'c4=']:
            for line in lines:
                if tag in line:
                    corners.append(map(float, line.split()[1:3]))

    # close the box by going back to origin
    corners.append(corners[0])
    # until each side has at least wanted number of points
    while len(corners) < 4 * min_edge_points:
        # work backwards, insertions don't change later indexes
        for i in xrange(len(corners) - 1, 0, -1):
            val = ll_mid(corners[i][0], corners[i][1], \
                    corners[i - 1][0], corners[i - 1][1])
            corners.insert(i, val)

    # write points the make the path
    with open(output, 'w') as mp:
        for point in corners:
            mp.write('%s %s\n' % (point[0], point[1]))


if __name__ == '__main__':
    # simple ability to process command line arguements
    import sys

    # option 1: give corners and output
    if len(sys.argv) >= 10:
        corners = []
        lonlats = map(float, sys.argv[1:9])
        for i in xrange(0, 8, 2):
            corners.append(lonlats[i:i + 2])
        output = sys.argv[9]

        # option 2: also specify min_edge_points
        if len(sys.argv) >= 11:
            points_per_edge = int(sys.argv[10])
            path_from_corners(corners = corners, output = output, \
                    min_edge_points = min_edge_points)
        else:
            path_from_corners(corners = corners, output = output)

    # option 3: give nothing, use velocity model and defaults
    else:
        path_from_corners()

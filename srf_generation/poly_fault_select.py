#!/usr/bin/env python
"""
Select faults based on being within a polygon.
"""

from argparse import ArgumentParser
import os

from matplotlib import pyplot as plt
import numpy as np
from osgeo import ogr, osr

from qcore.nhm import load_nhm


def points2polygon(points):
    """
    GDAL polygon from points (2D array).
    https://pcjericks.github.io/py-gdalogr-cookbook/geometry.html
    """
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for ll in np.atleast_2d(points):
        ring.AddPoint_2D(*map(float, ll))
    ring.CloseRings()
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    return poly


def points2trace(points):
    """
    GDAL trace from points (2D array).
    """
    trace = ogr.Geometry(ogr.wkbLineString)
    for ll in np.atleast_2d(points):
        trace.AddPoint_2D(*map(float, ll))
    return trace


if __name__ == "__main__":
    # arguments
    parser = ArgumentParser()
    parser.add_argument("nhm_file", help="path to NHM file", type=os.path.abspath)
    parser.add_argument(
        "points",
        help="longitude latitude of each point in polygon",
        nargs="+",
        type=float,
    )
    args = parser.parse_args()
    # process
    nhm = load_nhm(args.nhm_file)
    outline = np.array(args.points).reshape(-1, 2)

    # plot outline and store as GDAL object
    plt.plot(
        np.insert(outline[:, 0], 0, outline[-1, 0]),
        np.insert(outline[:, 1], 0, outline[-1, 1]),
        color="blue",
    )
    outline = points2polygon(outline)

    for fault in nhm.values():
        # store as GDAL object
        trace = points2trace(fault.trace)
        # intersection either point if point in linestring exact or linestring
        intersection = outline.Intersection(trace)
        if intersection.GetPointCount() == 0:
            # not intersecting
            colour = "black"
        else:
            # intersecting
            print(fault.name)
            if trace.Within(outline):
                colour = "red"
            else:
                colour = "orange"
        plt.plot(fault.trace[:, 0], fault.trace[:, 1], color=colour)
    plt.show()

#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example
python rel2vm_params.py ~/Data/Sources/Hossack/Srf/Hossack_REL01.csv ~/Data/VMs/Hossack
(if no output directory is specified, the same directory as REL01.csv will be used)

HELP:
python ./rel2vm_params.py -h
"""

from argparse import ArgumentParser
from glob import glob
from h5py import File as h5open
from logging import Logger
import math
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import platform
import subprocess

from shutil import rmtree, move

import sys
from tempfile import mkdtemp

from common import store_summary, temp_paths, NZ_LAND_OUTLINE, NZ_CENTRE_LINE

from srf_generation.input_file_generation.realisation_to_srf import get_corners_dbottom
from qcore import constants, geo, gmt, qclogging
from qcore.geo import R_EARTH
from qcore.utils import dump_yaml

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

from plot_vm import plot_vm

S_WAVE_KM_PER_S = 3.5

siteprop = Site()
siteprop.vs30 = 500
faultprop = Fault()
# DON'T CHANGE THIS - this assumes we have a surface point source
faultprop.ztor = 0.0

SPACE_LAND = 5.0 #min space between VM endge and land (km)
SPACE_SRF = 15.0 #min space between VM edge and SRF (km)
MIN_RJB = 0 # Specify a minimum horizontal distance (in km) for the VM to span
            # from the fault - invalid VMs will still not be generated",


# default scaling relationship
def mag2pgv(mag):
    return np.interp(
        mag,
        [3.5, 4.1, 4.7, 5.2, 5.5, 5.8, 6.2, 6.5, 6.8, 7.0, 7.4, 7.7, 8.0],
        [0.015, 0.0375, 0.075, 0.15, 0.25, 0.4, 0.7, 1.0, 1.35, 1.65, 2.1, 2.5, 3.0],
    )


# rrup at which pgv is close to target
def find_rrup(pgv_target):
    rrup = 50.0
    while True:
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        pgv = compute_gmm(faultprop, siteprop, GMM.Br_10, "PGV")[0]
        # factor 0.02 is conservative step to avoid infinite looping
        if pgv_target / pgv - 1 > 0.01:
            rrup -= rrup * 0.02
        elif pgv_target / pgv - 1 < -0.01:
            rrup += rrup * 0.02
        else:
            break
    return rrup, pgv


def determine_vm_extent(distance, hh, points, rot=0, wd="."):
    """
    rrup: in km
    """
    lon_mid, lat_mid, dx_km, dy_km = gmt.region_fit_oblique(points, 90 - rot, wd=wd)
    # extend by wanted rrup
    min_x = geo.ll_shift(lat_mid, lon_mid, distance + dx_km, 270 - rot)[::-1]
    max_x = geo.ll_shift(lat_mid, lon_mid, distance + dx_km, 90 - rot)[::-1]
    min_y = geo.ll_shift(lat_mid, lon_mid, distance + dy_km, 180 - rot)[::-1]
    max_y = geo.ll_shift(lat_mid, lon_mid, distance + dy_km, 0 - rot)[::-1]
    # mid, x, y extents
    return (
        (lon_mid, lat_mid),
        math.ceil(geo.ll_dist(min_x[0], min_x[1], max_x[0], max_x[1]) / hh) * hh,
        math.ceil(geo.ll_dist(min_y[0], min_y[1], max_y[0], max_y[1]) / hh) * hh,
    )


# return region that just fits the path made from 4 points
def corners2region(c1, c2, c3, c4):
    perimiter = np.array(geo.path_from_corners([c1, c2, c3, c4], output=None))
    x_min, y_min = np.min(perimiter, axis=0)
    x_max, y_max = np.max(perimiter, axis=0)
    return (x_min, x_max, y_min, y_max)


def build_corners(origin, rot, xlen, ylen):
    # wanted xlen, ylen is at corners
    # amount to shift from middle
    x_shift = xlen / 2.0
    y_shift = ylen / 2.0

    y_len_mid_shift = R_EARTH * math.asin(
        math.sin(y_shift / R_EARTH) / math.cos(x_shift / R_EARTH)
    )

    top_mid = geo.ll_shift(
        lat=origin[1], lon=origin[0], distance=y_len_mid_shift, bearing=rot
    )[::-1]
    bottom_mid = geo.ll_shift(
        lat=origin[1],
        lon=origin[0],
        distance=y_len_mid_shift,
        bearing=(rot + 180) % 360,
    )[::-1]

    top_mid_bearing = geo.ll_bearing(*top_mid, *origin) + 180 % 360
    bottom_mid_bearing = geo.ll_bearing(*bottom_mid, *origin)

    c2 = geo.ll_shift(
        lat=top_mid[1],
        lon=top_mid[0],
        distance=x_shift,
        bearing=(top_mid_bearing - 90) % 360,
    )[::-1]
    c1 = geo.ll_shift(
        lat=top_mid[1],
        lon=top_mid[0],
        distance=x_shift,
        bearing=(top_mid_bearing + 90) % 360,
    )[::-1]

    c3 = geo.ll_shift(
        lat=bottom_mid[1],
        lon=bottom_mid[0],
        distance=x_shift,
        bearing=(bottom_mid_bearing - 90) % 360,
    )[::-1]
    c4 = geo.ll_shift(
        lat=bottom_mid[1],
        lon=bottom_mid[0],
        distance=x_shift,
        bearing=(bottom_mid_bearing + 90) % 360,
    )[::-1]

    # at this point we have a perfect square (by corner distance)
    # c1 -> c4 == c2 -> c3 (right == left), c1 -> c2 == c3 -> c4 (top == bottom)
    return c1, c2, c3, c4


# proportion of gmt mask file filled (1s vs 0s)
# TODO: should be part of GMT library
def grd_proportion(grd_file):
    try:
        with h5open(grd_file, "r") as grd:
            return 100.0 * np.nansum(grd["z"][...]) / float(grd["z"].size)
    except IOError:
        print("INVALID GRD: %s" % (grd_file))
        # this should no longer happen when auto_dxy is used
        raise


# keep dx and dy small enough relative to domain
# h5py will crash when opening grd with less than 2^14 gridpoints
def auto_dxy(region):
    # assume not going over earth quadrants
    x_diff = region[1] - region[0]
    y_diff = region[3] - region[2]

    # start off at roughly 1km
    dxy = 0.01
    # this formula may not be exact, lower value than actual if anything (safe)
    while math.ceil(x_diff / dxy) * math.ceil(y_diff / dxy) < 16384:
        dxy *= 0.5
    return dxy


# z extent based on magnitude, hypocentre depth
def auto_z(mag, depth):
    return round(
        10
        + depth
        + (10 * np.power((0.5 * np.power(10, (0.55 * mag - 1.2)) / depth), (0.3))),
        0,
    )


# simulation time based on area
def auto_time2(
    vm_corners,
    srf_corners,
    ds_multiplier,
    depth=0,
    logger: Logger = qclogging.get_basic_logger(),
):
    """Calculates the sim duration from the bounds of the vm and srf
    :param vm_corners: A numpy array of (lon, lat) tuples
    :param srf_corners: A [4n*2] numpy array of (lon, lat) pairs where n is the number of planes in the srf
    :param ds_multiplier: An integer
    :param logger: The logger to pass all messages to
    :param depth: Depth of the fault. This allows the horizontal distance to be extended to include the distance to
    travel to the surface too"""
    # S wave arrival time is determined by the distance from the srf centroid to the furthest corner
    h_dist = geo.get_distances(
        vm_corners,
        (srf_corners[:, 0].max() + srf_corners[:, 0].min()) / 2,
        (srf_corners[:, 1].max() + srf_corners[:, 1].min()) / 2,
    ).max()
    s_wave_arrival = (h_dist ** 2 + depth ** 2) ** 0.5 / S_WAVE_KM_PER_S
    logger.debug("s_wave_arrival: {}".format(s_wave_arrival))
    # Rrup is determined by the largest vm corner to nearest srf corner distance
    rjb = max([geo.get_distances(srf_corners, *corner).min() for corner in vm_corners])
    siteprop.Rrup = (rjb ** 2 + depth ** 2) ** 0.5
    logger.debug("rrup: {}".format(siteprop.Rrup))
    # magnitude is in faultprop
    ds = compute_gmm(faultprop, siteprop, GMM.AS_16, "Ds595")[0]
    logger.debug("ds: {}".format(ds))
    logger.debug("ds_multiplier: {}".format(ds_multiplier))
    return s_wave_arrival + ds_multiplier * ds


def vm_land(c1, c2, c3, c4, wd="."):
    """
    Proportion of VM on land.
    """
    path_vm = "%s/mask_vm.path" % (wd)
    grd_vm = "%s/mask_vm.grd" % (wd)
    grd_land = "%s/mask_land.grd" % (wd)
    grd_vm_land = "%s/mask_vm_land.grd" % (wd)

    geo.path_from_corners([c1, c2, c3, c4], output=path_vm)
    vm_region = corners2region(c1, c2, c3, c4)
    dxy = auto_dxy(vm_region)

    gmt.grd_mask("f", grd_land, region=vm_region, dx=dxy, dy=dxy, wd=wd)
    gmt.grd_mask(path_vm, grd_vm, region=vm_region, dx=dxy, dy=dxy, wd=wd)
    gmt.grdmath([grd_land, grd_vm, "BITAND", "=", grd_vm_land], wd=wd)

    return grd_proportion(grd_vm_land) / grd_proportion(grd_vm) * 100


# get longitude on NZ centre path at given latitude
def centre_lon(lat_target):
    # find closest points either side
    # TODO: use np.interp
    with open(NZ_CENTRE_LINE, "r") as c:
        for line in c:
            ll = list(map(float, line.split()))
            if ll[1] < lat_target:
                ll_bottom = ll
            else:
                ll_top = ll
                break
    try:
        # find rough proportion of distance
        lat_ratio = (lat_target - ll_bottom[1]) / (ll_top[1] - ll_bottom[1])
    except (UnboundLocalError, NameError):
        print("ERROR: VM edge out of lat bounds for centre line.")
        raise
    # consider same ratio along longitude difference to be correct
    return ll_bottom[0] + lat_ratio * (ll_top[0] - ll_bottom[0])


# maximum width along lines a and b
def reduce_domain(
    origin,
    bearing,
    xlen,
    ylen,
    hh,
    space_srf = SPACE_SRF,
    space_land = SPACE_LAND,
    wd,
    logger: Logger = qclogging.get_basic_logger(),
):
    """Reduces the domain of the VM by removing areas that are over sea only. Gives a buffer around coast line and the
    srf boundary. Returns the new values to be used.
    Uses great circle geometry"""
    # number of scan lines accross domain (inclusive of edges)
    scanlines = max(round(ylen / 5), 81)

    # store scan data for 2nd level processing
    scan_extremes = np.zeros((scanlines, 2, 2))

    # Half the x/ylen, commonly used value
    x_shift = xlen / 2.0
    y_shift = ylen / 2.0

    # The first point will always change them
    # Positive is east of mid line
    dist_east_from_mid = -x_shift
    # Positive is west of mid line
    dist_west_from_mid = -x_shift

    # The distance between the origin and domain boundary mid points
    y_len_mid_shift = R_EARTH * math.asin(
        math.sin(y_shift / R_EARTH) / math.cos(x_shift / R_EARTH)
    )
    x_len_mid_shift = R_EARTH * math.asin(
        math.sin(x_shift / R_EARTH) / math.cos(y_shift / R_EARTH)
    )

    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint=True)):
        m = geo.ll_shift(
            origin[1], origin[0], y_len_mid_shift * 2 * (x - 0.5), bearing
        )[::-1]
        # extend path a -> b, make sure we reach land for intersections
        bearing_p = geo.ll_bearing(m[0], m[1], origin[0], origin[1])
        if 90 <= bearing_p <= 270:
            # Make sure bearing is pointing northwards: [-90, 90] % 360
            bearing_p = (bearing_p + 180) % 360
        ap = geo.ll_shift(m[1], m[0], 600 + x_len_mid_shift, (bearing_p + 90) % 360)[
            ::-1
        ]
        bp = geo.ll_shift(m[1], m[0], 600 + x_len_mid_shift, (bearing_p + 270) % 360)[
            ::-1
        ]

        # GMT spatial is not great-circle path connective
        geo.path_from_corners(
            corners=[ap, bp],
            output="%s/tempEXT.tmp" % (wd),
            min_edge_points=80,
            close=False,
        )
        isections, icomps = gmt.intersections(
            ["%s/srf.path" % (wd), NZ_LAND_OUTLINE, "%s/tempEXT.tmp" % (wd)],
            items=True,
            containing="%s/tempEXT.tmp" % (wd),
        )
        if len(isections) == 0:
            scan_extremes[i] = np.nan
            continue
        # shift different intersections by item-specific padding amount

        for c, intersection in enumerate(isections):
            if "%s/srf.path" % (wd) in icomps[c]:
                diff = space_srf
            elif NZ_LAND_OUTLINE in icomps[c]:
                diff = space_land
            else:
                # Something went wrong
                continue

            isection_bearing = geo.ll_bearing(*m, *intersection)
            if (
                abs(isection_bearing - (bearing_p + 90)) < 5
                and dist_east_from_mid < x_shift
            ):
                # The intersection is (relatively) east of the mid point
                east = geo.ll_dist(*intersection, *m) + diff
                west = diff - geo.ll_dist(*intersection, *m)
            elif (
                abs(isection_bearing - (bearing_p + 270) % 360) < 5
                and dist_west_from_mid < x_shift
            ):
                # The intersection is (relatively) west of the mid point
                east = diff - geo.ll_dist(*intersection, *m)
                west = geo.ll_dist(*intersection, *m) + diff
            else:
                # If the point isn't east or west, it's probably on the center line
                east = west = diff

            # Only update if we have a better choice, but only go up to the original value
            if east > dist_east_from_mid:
                dist_east_from_mid = min(east, x_shift)
            if west > dist_west_from_mid:
                dist_west_from_mid = min(west, x_shift)

        # for north/south (later on), use direct/un-shifted points
        scan_extremes[i] = (
            isections[np.argmax(isections, axis=0)[0]],
            isections[np.argmin(isections, axis=0)[0]],
        )

    # second scan, reduce north and south
    over_s = 0
    over_n = ylen
    land_discovered = False
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint=True)):
        m = geo.ll_shift(
            origin[1], origin[0], y_len_mid_shift * 2 * (x - 0.5), bearing
        )[::-1]
        # extend path a -> b, make sure we reach land for intersections
        bearing_p = geo.ll_bearing(m[0], m[1], origin[0], origin[1])
        if 90 <= bearing_p <= 270:
            bearing_p = (bearing_p + 180) % 360
        # Extend scan points to the relative east and west, at the distance of the widest part.
        # TODO: use the width of the domain at that point instead of the max width
        a = geo.ll_shift(m[1], m[0], x_len_mid_shift, bearing_p + 90)[::-1]
        b = geo.ll_shift(m[1], m[0], x_len_mid_shift, bearing_p + 270)[::-1]
        # determine if window contains any land by locations of intersections
        isect_inside = False
        isect_east = False
        isect_west = False
        for extreme in scan_extremes[i]:
            if b[0] <= extreme[0] <= a[0]:
                isect_inside = True
            elif extreme[0] > a[0]:
                isect_east = True
            else:
                isect_west = True
        # not on land if isect outside zone and either no east or no west isect
        if isect_inside or (isect_east and isect_west):
            if not land_discovered:
                # shift bottom edge up
                over_s = math.floor((x * ylen - 15) / hh) * hh
            land_discovered = True
            over_n = math.floor((ylen - x * ylen - 15) / hh) * hh

    if not np.isclose(dist_east_from_mid, dist_west_from_mid):
        # If the east and west are not changing by the same amount, then we have to move the origin, and get the new
        # bearing
        logger.debug(
            "Pruning more from East or West than the other, moving center point"
        )
        x_mid_ratio = x_len_mid_shift / x_shift
        origin1 = geo.ll_shift(
            *origin[::-1],
            x_mid_ratio * (dist_east_from_mid - dist_west_from_mid) / 2,
            bearing + 90,
        )[::-1]
        bearing = geo.ll_bearing(*origin1, *origin)
        if dist_east_from_mid > dist_west_from_mid:
            # If we moved east we need the right angle clockwise from the old origin
            bearing = (bearing + 90) % 360
        else:
            # If we moved west then we need the right angle anticlockwise from the old origin
            bearing = (bearing - 90) % 360
        origin = origin1
        logger.debug(
            "New origin after x shift: {}, new bearing after x shift: {}".format(
                origin, bearing
            )
        )

    xlen = math.ceil((dist_east_from_mid + dist_west_from_mid) / hh) * hh
    logger.debug("New xlen: {}".format(xlen))

    # If either of the numbers is negative then we don't trim that edge, so set the minimum value of each to 0
    over_n_max = max(over_n, 0)
    over_s_max = max(over_s, 0)
    if not np.isclose(over_n_max, over_s_max):
        logger.debug(
            "Pruning more from North or South than the other, moving center point"
        )
        # If we don't move the north and south boundaries in by the same amount, then we need to move the origin an
        # amount relative to the change
        # The multiplier for the difference between the origin->mid north point and the mid east -> northeast corner
        # This is symmetric with north/south and east/west
        y_mid_ratio = y_len_mid_shift / y_shift
        # Move the origin by a distance relative to the amount taken off each end, multiplied by the ratio
        origin1 = geo.ll_shift(
            *origin[::-1], y_mid_ratio * (over_s_max - over_n_max) / 2, bearing
        )[::-1]
        bearing = geo.ll_bearing(*origin1, *origin)
        if over_s_max > over_n_max:
            # The point moved north, so we are looking South
            bearing = (bearing + 180) % 360
        origin = origin1
        logger.debug(
            "New origin after y shift: {}, new bearing after y shift: {}".format(
                origin, bearing
            )
        )

    ylen -= over_s_max + over_n_max
    logger.debug("New ylen: {}".format(ylen))

    return origin, bearing, xlen, ylen


def optimise_vm_params(
    srf_meta, options, out_dir, ptemp, logger: Logger = qclogging.get_basic_logger()
):
    """
    Optimises VM domain and returns a dictionary that extends vm_params_dict

    Parameters
    ----------
    srf_meta : Data extracted from REL csv file
    options : Default/User-specified options (via arguments)
    out_dir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    ptemp : Temporary work directory (eg. .../Data/VMs/Hossack/_tmp_Hossack_nho8ui41)
    logger :

    Returns
    -------
    vm_params_dict_extended that contains extra info not stored in vm_params.yaml

    """
    # properties stored in classes (fault of external code)
    faultprop.Mw = srf_meta["mag"]
    faultprop.rake = srf_meta["rake"]
    faultprop.dip = srf_meta["dip"]
    faultprop.faultstyle = None
    # rrup to reach wanted PGV
    if options["pgv"] == -1.0:
        pgv = mag2pgv(faultprop.Mw)
        logger.debug("pgv determined from magnitude. pgv set to {}".format(pgv))
    else:
        pgv = options["pgv"]
        logger.debug("pgv set to {}".format(pgv))
    rrup, pgv_actual = find_rrup(pgv)

    if "dtop" in srf_meta:
        fault_depth = np.min(srf_meta["dtop"])
    else:
        fault_depth = srf_meta["hdepth"]

    rjb = 0
    if fault_depth < rrup * 2:
        # rjb = (rrup ** 2 - fault_depth ** 2) ** 0.5
        rjb = max(
            MIN_RJB, rrup
        )  # sets rrup equal to rjb to ensure deep ruptures have sufficient VM size

    # original, unrotated vm
    bearing = 0
    origin, xlen0, ylen0 = determine_vm_extent(
        rjb, options["hh"], srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=ptemp
    )

    if fault_depth > rrup:
        logger.warning(
            "fault_depth > rrup. Fault too deep, will not generate VM. "
            f"Setting xlen, ylen to 0. Previous values: x:{xlen0}, y:{ylen0}"
        )
        xlen0 = 0
        ylen0 = 0

    o1, o2, o3, o4 = build_corners(origin, bearing, xlen0, ylen0)
    vm0_region = corners2region(o1, o2, o3, o4)
    plot_region = (
        vm0_region[0] - 1,
        vm0_region[1] + 1,
        vm0_region[2] - 1,
        vm0_region[3] + 1,
    )

    # proportion in ocean
    if options["no_optimise"] or xlen0 <= 0 or ylen0 <= 0:
        logger.debug("Not optimising for land coverage")
        land0 = 100
    else:
        land0 = vm_land(o1, o2, o3, o4, wd=ptemp)
        logger.debug("Land coverage found to be {}%".format(land0))

    # for plotting and calculating VM domain distance
    with open("%s/srf.path" % (ptemp), "wb") as sp:
        for plane in srf_meta["corners"]:
            sp.write("> srf plane\n".encode())
            np.savetxt(sp, plane, fmt="%f")
            sp.write("%f %f\n".encode() % (tuple(plane[0])))

    # modify VM if necessary
    adjusted = False
    if faultprop.Mw >= 3.5 and land0 < 99:
        adjusted = True
        logger.info(
            "Fault {} has magnitude greater than 3.5 and less than 99% land coverage. Optimising.".format(
                srf_meta["name"]
            )
        )

        # rotation based on centre line bearing at this point
        l1 = centre_lon(vm0_region[2])
        l2 = centre_lon(vm0_region[3])
        mid = geo.ll_mid(l1, vm0_region[2], l2, vm0_region[3])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, vm0_region[3]))

        # wanted distance is at corners, not middle top to bottom
        _, xlen1, ylen1 = determine_vm_extent(
            rjb,
            options["hh"],
            srf_meta["corners"].reshape((-1, 2)),
            rot=bearing,
            wd=ptemp,
        )

        logger.debug(
            "Pre-optimisation origin: {}, bearing: {}, xlen: {}, ylen: {}".format(
                origin, bearing, xlen1, ylen1
            )
        )
        # cut down ocean areas
        origin, bearing, xlen1, ylen1, = reduce_domain(
            origin,
            bearing,
            xlen1,
            ylen1,
            options["hh"],
            ptemp,
            logger=logger,
        )
        logger.debug(
            "After optimisation origin: {}, bearing: {}, xlen: {}, ylen: {}".format(
                origin, bearing, xlen1, ylen1
            )
        )

        try:
            c1, c2, c3, c4 = build_corners(origin, bearing, xlen1, ylen1)
        except ValueError:
            error_message = "Error for vm {}. Raising exception.".format(
                srf_meta["name"]
            )
            logger.log(qclogging.NOPRINTCRITICAL, error_message)
            raise ValueError(error_message)
        else:
            logger.debug("Optimised corners of the VM are: {}".format((c1, c2, c3, c4)))

        # proportion in ocean
        land1 = vm_land(c1, c2, c3, c4, wd=ptemp)
        logger.debug("Optimised land coverage found to be {}%".format(land1))

        # adjust region to fit new corners
        plot_region = (
            min(c4[0], c1[0], c3[0], c2[0], plot_region[0]),
            max(c4[0], c1[0], c3[0], c2[0], plot_region[1]),
            min(c4[1], c1[1], c3[1], c2[1], plot_region[2]),
            max(c4[1], c1[1], c3[1], c2[1], plot_region[3]),
        )
    else:
        logger.info(
            "Magnitude less than 3.5 or land coverage is greater than 99%, not modifying {}".format(
                srf_meta["name"]
            )
        )
        # store original variables as final versions
        ylen1 = ylen0
        xlen1 = xlen0
        land1 = land0
        c1, c2, c3, c4 = o1, o2, o3, o4
    # not enough land in final domain
    if math.floor(land1) == 0:
        logger.info(
            "Land coverage is less than 1%. Setting xlen and ylen to 0, not creating VM"
        )
        xlen1 = 0
        ylen1 = 0

    # zlen is independent from xlen and ylen
    zlen = (
        round(auto_z(faultprop.Mw, srf_meta["dbottom"]) / options["hh"]) * options["hh"]
    )
    logger.debug("zlen set to {}".format(zlen))

    if xlen1 != 0 and ylen1 != 0 and zlen != 0:

        # modified sim time
        vm_corners = np.asarray([c1, c2, c3, c4])
        initial_time = auto_time2(
            vm_corners,
            np.concatenate(srf_meta["corners"], axis=0),
            options["ds_multiplier"],
            fault_depth,
            logger=logger,
        )
        sim_time1 = (initial_time // options["dt"]) * options["dt"]
        logger.debug(
            "Unrounded sim duration: {}. Rounded sim duration: {}".format(
                initial_time, sim_time1
            )
        )

        # optimisation results
        return {
            "name": srf_meta["name"],
            "mag": faultprop.Mw,
            "dbottom": srf_meta["dbottom"],
            "zlen": zlen,
            "sim_time": sim_time1,
            "xlen": xlen0,
            "ylen": ylen0,
            "land": land0,
            "zlen_mod": zlen,
            "sim_time_mod": sim_time1,
            "xlen_mod": xlen1,
            "ylen_mod": ylen1,
            "land_mod": land1,
            "origin": origin,
            "bearing": bearing,
            "adjusted": adjusted,
            "plot_region": plot_region,
            "path": "%s %s\n%s %s\n%s %s\n%s %s\n"
            % (o1[0], o1[1], o2[0], o2[1], o3[0], o3[1], o4[0], o4[1]),
            "path_mod": "%s %s\n%s %s\n%s %s\n%s %s\n"
            % (c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]),
        }

    return None  # failed to create VM if it is entirely in the ocean


def save_vm_params_yaml(
    vm_params_path=None,
    vm_dir=None,
    origin=(170, -40),
    rot=0,
    xlen=100,
    ylen=100,
    zmax=40,
    zmin=0,
    hh=0.4,
    min_vs=0.5,
    mag=5.5,
    centroid_depth=7,
    sim_duration=100,
    code="rt",
    model_version="1.65",
    topo_type="BULLDOZED",
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Save vm_params.yaml with the specified info

    Parameters
    ----------
    vm_params_path : Path to vm_params.yaml
    vm_dir : A relative path (default: output) where NZVM will place output files (cannot exist)
    origin :
    rot :
    xlen :
    ylen :
    zmax :
    zmin :
    hh :
    min_vs :
    mag :
    centroid_depth :
    sim_duration :
    code :
    model_version :
    topo_type :
    logger :
    """

    if vm_params_path is not None:
        # must also convert mag from np.float to float
        vm_params_dict = {
            "mag": float(mag),
            "centroidDepth": float(centroid_depth),
            "MODEL_LAT": float(origin[1]),
            "MODEL_LON": float(origin[0]),
            "MODEL_ROT": float(rot),
            "hh": hh,
            "min_vs": float(min_vs),
            "model_version": model_version,
            "topo_type": topo_type,
            "output_directory": os.path.basename(vm_dir),
            "extracted_slice_parameters_directory": "SliceParametersNZ/SliceParametersExtracted.txt",
            "code": code,
            "extent_x": float(xlen),
            "extent_y": float(ylen),
            "extent_zmax": float(zmax),
            "extent_zmin": float(zmin),
            "sim_duration": float(sim_duration),
            "flo": min_vs / (5.0 * hh),
            "nx": int(round(float(xlen) / hh)),
            "ny": int(round(float(ylen) / hh)),
            "nz": int(round(float(zmax - zmin) / hh)),
            "sufx": "_%s01-h%.3f" % (code, hh),
        }

        dump_yaml(
            vm_params_dict,
            "{}".format(vm_params_path),
        )
        logger.debug("Saved {}".format(vm_params_path))
        return vm_params_dict
    return None


def main(
    srf_meta, options, out_dir, logger_name: str = "rel2vm_params", plot_enabled=True
):
    """
    Orchestrates conversion from rel CSV file to vm_params.yaml
    Parameters
    ----------
    srf_meta : Data extracted from REL csv file
    options : Default/User-specified options (via arguments)
    out_dir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    logger_name :
    plot_enabled :

    Returns
    -------
    vm_params_dict: Dictionary already dumped to vm_params.yaml, Returned for store_summary
    """

    logger = qclogging.get_realisation_logger(
        qclogging.get_logger(logger_name), srf_meta["name"]
    )

    
    # temp directory for current process
    ptemp = mkdtemp(prefix="_tmp_%s_" % (srf_meta["name"]), dir=out_dir)

    vm_params_dict = {}
    vm_params_dict_extended = optimise_vm_params(
        srf_meta, options, out_dir, ptemp, logger=logger
    )
    if vm_params_dict_extended is not None:
        # store configs
        vm_working_dir, _, vm_params_path = temp_paths(ptemp)

        vm_params_dict = save_vm_params_yaml(
            vm_params_path=vm_params_path,
            vm_dir=vm_working_dir,
            origin=vm_params_dict_extended["origin"],
            rot=vm_params_dict_extended["bearing"],
            xlen=vm_params_dict_extended["xlen_mod"],
            ylen=vm_params_dict_extended["ylen_mod"],
            zmax=vm_params_dict_extended["zlen_mod"],
            hh=options["hh"],
            min_vs=options["min_vs"],
            mag=faultprop.Mw,
            centroid_depth=srf_meta["hdepth"],
            sim_duration=vm_params_dict_extended["sim_time_mod"],
            topo_type=options["vm_topo"],
            model_version=options["vm_version"],
        )
        if vm_params_dict:
            # save important files
            logger.debug("Creating directory {}".format(vm_working_dir))
            os.makedirs(vm_working_dir)
            # move(nzvm_cfg_path, vm_working_dir)
            vm_params_path_final = os.path.join(
                out_dir, os.path.basename(vm_params_path)
            )
            if os.path.exists(vm_params_path_final):
                logger.info(
                    f"Warning: {vm_params_path_final} already exists and to be overwritten."
                )
                os.remove(vm_params_path_final)

            move("{}".format(vm_params_path), out_dir)
            logger.info("Wrote vm_params.yaml at {}".format(out_dir))
            if plot_enabled:
                plot_vm(
                    vm_params_dict_extended,
                    srf_meta["corners"],
                    faultprop.Mw,
                    out_dir,
                    ptemp,
                    logger=logger,
                )

            print(f"Success: Wrote vm_params.yaml at {out_dir}", file=sys.stderr)
            # generate a corners like NZVM would have
            logger.debug("Saving VeloModCorners.txt")
            with open("{}/VeloModCorners.txt".format(out_dir), "wb") as c:
                c.write("> VM corners (python generated)\n".encode())
                c.write(">Lon    Lat\n".encode())
                c.write(vm_params_dict_extended["path_mod"].encode())

    # working dir cleanup, return info about VM
    logger.debug("Cleaning up temp directory {}".format(ptemp))
    rmtree(ptemp)
    return vm_params_dict


def load_msgs(
        ds_multiplier,
        dt,
        hh,
        min_vs,
        no_optimise,
        out_dir,
        pgv,
        rel_file,
        vm_topo,
        vm_version,
        logger: Logger = qclogging.get_basic_logger()
):
    # returns list of appropriate srf metadata
    msgs = []

    #rel_file is assumed to be the first realisation
    if not os.path.exists(rel_file):
        logger.info("REL csv file not found: {}".format(rel_file))
        return

    # name is unique and based on basename
    name = os.path.splitext(os.path.basename(rel_file))[0].split("_")[0]
    logger.debug("Found first REL file for {}".format(name))

    a = pd.read_csv(rel_file)
    rake = a["rake"].loc[0]

    planes = []
    num_subfaults = a["plane_count"].loc[0]
    for i in range(num_subfaults):
        plane = {}
        plane["centre"] = [
            a["clon_subfault_{}".format(i)].loc[0],
            a["clat_subfault_{}".format(i)].loc[0],
        ]
        plane["length"] = a["length_subfault_{}".format(i)].loc[0]
        plane["width"] = a["width_subfault_{}".format(i)].loc[0]
        plane["strike"] = a["strike_subfault_{}".format(i)].loc[0]
        plane["dip"] = a["dip_subfault_{}".format(i)].loc[0]
        plane["dtop"] = a["dtop_subfault_{}".format(i)].loc[0]
        planes.append(plane)

    corners, dbottom = get_corners_dbottom(planes, dip_dir=a["dip_dir"].loc[0])

    mag = a["magnitude"].loc[0]
    # try:
    #    mag = mom2mag(sum(map(mag2mom, a["magnitude"].loc[0])))
    # except TypeError:
    #    mag = a["magnitude"].loc[0]
    dtop = a["dtop"].loc[0]
    dbottom = a["dbottom"].loc[0]

    hdepth = (
        np.round(a["dhypo"].loc[0], decimals=1)
        * np.sin(np.radians(a["dip"].loc[0]))
        + dtop
    )

    msgs.append(
        (
            {
                "name": name,
                "dip": a["dip"].loc[0],
                "rake": rake,
                "dbottom": dbottom,
                "corners": np.array(corners),
                "mag": mag,
                "hdepth": hdepth,
            },
            {
                "pgv": pgv,
                "hh": hh,
                "no_optimise": no_optimise,
                "ds_multiplier": ds_multiplier,
                "dt": dt,
                "min_vs": min_vs,
                "vm_topo": vm_topo,
                "vm_version": vm_version,
            },
            out_dir,
            qclogging.get_realisation_logger(logger, name).name,
            True,
        )
    )

    return msgs


def load_args(logger: Logger = qclogging.get_basic_logger()):
    parser = ArgumentParser()
    arg = parser.add_argument

    arg("rel_file",help="REL csv file")

    arg("--out_dir", help="output directory to place VM files "
                          "(if not specified, the same path as rel_file is used", default=None)

    arg(
        "--pgv",
        help="max PGV at velocity model perimiter (estimated, cm/s)",
        type=float,
        default=5.0,
    )
    arg("--hh", help="velocity model grid spacing (km)", type=float, default=0.4)
    arg(
        "--dt",
        help="timestep to estimate simulation duration (s) Default: hh/20",
        type=float,
        default=None,
    )

    arg("--min-vs", help="for nzvm gen and flo (km/s)", type=float, default=0.5)

    arg("--vm-version", help="velocity model version to generate", default="2.06")
    arg(
        "--vm-topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )
    arg(
        "--no_optimise",
        help="Don't try and optimise the vm if it is off shore. Removes dependency on having GMT coastline data",
        action="store_true",
    )
    arg(
        "--min-rjb",
        help="Specify a minimum horizontal distance (in km) for the VM to span from the fault"
        " - invalid VMs will still not be generated",
        default=0,
    )
    arg(
        "--ds-multiplier",
        help="Sets the DS multiplier for setting the sim-duration. Validation runs default to 1.2. Cybershake runs"
        "should manually set it to 0.75",
        default=1.2,
        type=float,
    )
    args = parser.parse_args()
    if args.dt is None:
        args.dt = args.hh / 20

    if args.out_dir is None:
        args.out_dir = os.path.abspath(os.path.dirname(args.rel_file))


    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("rel2vm_params")
    qclogging.add_general_file_handler(
        logger, os.path.join(os.getcwd(), "rel2vm_params_log.txt")
    )
    args = load_args(logger=logger)

    logger.debug("Loading REL csv")
    msg_list = load_msgs(args.ds_multiplier,args.dt,args.hh, args.min_vs,
                         args.no_optimise,args.out_dir,args.pgv,args.vm_topo,args.vm_version,
                         logger=logger)

    if len(msg_list) == 0:
        message = "Found nothing to do, exiting."
        logger.log(qclogging.NOPRINTCRITICAL, message)
        sys.exit(message)

    # prepare to run
    os.makedirs(args.out_dir, exist_ok=True)

    reports = main(*msg_list)

    # store summary
    store_summary(
        os.path.join(args.out_Dir, "rel2vm_params_info.csv"), reports, logger=logger
    )

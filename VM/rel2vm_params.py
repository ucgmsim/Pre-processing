#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example
python rel2vm_params.py ~/Data/Sources/Hossack/Srf/Hossack_REL01.csv -o ~/Data/VMs/Hossack
(if no output directory is specified, the same directory as REL01.csv will be used)

HELP:
python ./rel2vm_params.py -h
"""

from argparse import ArgumentParser
from h5py import File as h5open
from logging import Logger
import math
import numpy as np
import os
import pandas as pd
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from srf_generation.input_file_generation.realisation_to_srf import get_corners_dbottom

from qcore import constants, geo, gmt, qclogging
from qcore.geo import R_EARTH
from qcore.utils import dump_yaml
from qcore.simulation_structure import get_fault_from_realisation
from qcore.validate_vm import validate_vm_bounds

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

from plot_vm import plot_vm

script_dir = Path(__file__).resolve().parent
NZ_CENTRE_LINE = script_dir / "../SrfGen/NHM/res/centre.txt"
NZ_LAND_OUTLINE = script_dir / "../SrfGen/NHM/res/rough_land.txt"

S_WAVE_KM_PER_S = 3.5

siteprop = Site()
siteprop.vs30 = 500
faultprop = Fault()
# DON'T CHANGE THIS - this assumes we have a surface point source
faultprop.ztor = 0.0

SPACE_LAND = 5.0  # min space between VM edge and land (km)
SPACE_SRF = 15.0  # min space between VM edge and SRF (km)
MIN_RJB = (
    0
)  # minimum horizontal distance (in km) for the VM to span from the fault - invalid VMs will still not be generated",


#
def mag2pgv(mag: float):
    """
    Return PGV for a given magnitude based on default scaling relationship
    Parameters
    ----------
    mag : magnitude

    Returns
    -------
    Interpolated pgv

    """
    return np.interp(
        mag,
        [3.5, 4.1, 4.7, 5.2, 5.5, 5.8, 6.2, 6.5, 6.8, 7.0, 7.4, 7.7, 8.0],
        [0.015, 0.0375, 0.075, 0.15, 0.25, 0.4, 0.7, 1.0, 1.35, 1.65, 2.1, 2.5, 3.0],
    )


def find_rrup(pgv_target: float):
    """
    Find rrup at which pgv is close to target. Has to be computed iteratively as GMPE is irreversible.

    Parameters
    ----------
    pgv_target : pgv that we wish to come close to

    Returns
    -------

    """
    rrup = 50.0
    pgv = np.inf
    while pgv_target / pgv < 0.99 or pgv_target / pgv > 1.01:
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        pgv = compute_gmm(faultprop, siteprop, GMM.Br_10, "PGV")[0]
        # factor 0.02 is conservative step to avoid infinite looping
        if pgv_target / pgv > 1.01:
            rrup -= rrup * 0.02
        elif pgv_target / pgv < 0.99:
            rrup += rrup * 0.02

    return rrup, pgv


def determine_vm_extent(
    distance: float, hh: float, points: np.ndarray, rot: float = 0, wd: Path = Path(".")
):
    """
    Determine the extent of velocity model

    Parameters
    ----------
    distance :
    hh :
    points : array of points in (lon,lat) format
    rot :
    wd : working directory

    Returns
    -------
    ( (lan_mid,lat_mid), x extents, y extents )

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


def corners2region(c1: tuple, c2: tuple, c3: tuple, c4: tuple):
    """
    return region that just fits the path made from 4 points
    Parameters
    ----------
    c1 : tuple of (lon,lat)
    c2 : as above
    c3 : as above
    c4 : as above

    Returns
    -------
    (x_min, x_max, y_min, y_max) : (min lon, max lon, min lat, max lat)

    """
    perimeter = np.array(geo.path_from_corners([c1, c2, c3, c4], output=None))
    x_min, y_min = np.min(perimeter, axis=0)
    x_max, y_max = np.max(perimeter, axis=0)
    return (x_min, x_max, y_min, y_max)


def grd_proportion(grd_file: Path):
    """
    Get proportion of gmt mask file filled (1s vs 0s)

    TODO: should be part of GMT library
    Parameters
    ----------
    grd_file :

    Returns
    -------
    Proportion (%) of gmt mask file filled

    """
    try:
        with h5open(grd_file, "r") as grd:
            return 100.0 * np.nansum(grd["z"][...]) / float(grd["z"].size)
    except IOError:
        print(f"INVALID GRD: {grd_file}")
        # this should no longer happen when get_dxy is used
        raise


def get_dxy(region: tuple):
    """
    Get the grid spacing size that will be used to open a grd file
    Note: h5py will crash when opening grd with less than 2^14 gridpoints,
    so keep dx and dy small enough relative to domain

    Parameters
    ----------
    region : (min lon, max lon, min lat, max lat)

    Returns
    -------
    dxy

    """
    # assume not going over earth quadrants
    x_diff = region[1] - region[0]
    y_diff = region[3] - region[2]

    # start off at roughly 1km
    dxy = 0.01
    # this formula may not be exact, lower value than actual if anything (safe)
    while math.ceil(x_diff / dxy) * math.ceil(y_diff / dxy) < 16384:
        dxy *= 0.5
    return dxy


def get_max_depth(mag: float, depth: float):
    """
    Maximum depth (ie. z extent) based on magnitude, hypocentre depth
    Parameters
    ----------
    mag : magnitude
    depth : hypocentre depth

    Returns
    -------
    maximum depth
    """
    return round(
        10
        + depth
        + (10 * np.power((0.5 * np.power(10, (0.55 * mag - 1.2)) / depth), (0.3))),
        0,
    )


def get_sim_duration(
    vm_corners: np.ndarray,
    srf_corners: np.ndarray,
    ds_multiplier: int,
    depth: float = 0,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Calculates the sim duration from the bounds of the vm and srf

    Parameters
    ----------
    vm_corners : A numpy array of (lon, lat) tuples
    srf_corners : A [4n*2] numpy array of (lon, lat) pairs where n is the number of planes in the srf
    ds_multiplier :
    depth : Depth of the fault. This allows the horizontal distance to be extended to include the distance to
    travel to the surface too
    logger :

    Returns
    -------
    time in seconds
    """
    # S wave arrival time is determined by the distance from the srf centroid to the furthest corner
    h_dist = geo.get_distances(
        vm_corners,
        (srf_corners[:, 0].max() + srf_corners[:, 0].min()) / 2,
        (srf_corners[:, 1].max() + srf_corners[:, 1].min()) / 2,
    ).max()
    s_wave_arrival = (h_dist ** 2 + depth ** 2) ** 0.5 / S_WAVE_KM_PER_S
    logger.debug(f"s_wave_arrival: {s_wave_arrival}")
    # Rrup is determined by the largest vm corner to nearest srf corner distance
    rjb = max([geo.get_distances(srf_corners, *corner).min() for corner in vm_corners])
    siteprop.Rrup = (rjb ** 2 + depth ** 2) ** 0.5
    logger.debug(f"rrup: {siteprop.Rrup}")
    # magnitude is in faultprop
    ds = compute_gmm(faultprop, siteprop, GMM.AS_16, "Ds595")[0]
    logger.debug(f"ds: {ds}")
    logger.debug(f"ds_multiplier: {ds_multiplier}")
    return s_wave_arrival + ds_multiplier * ds


def reduce_domain(
    origin: float,
    bearing: float,
    xlen: float,
    ylen: float,
    hh: float,
    wd: Path,
    space_srf: float = SPACE_SRF,
    space_land: float = SPACE_LAND,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Reduces the domain of the VM by removing areas that are over sea only. Gives a buffer around coast line and the
    srf boundary.
    Uses great circle geometry

    Parameters
    ----------
    origin :
    bearing :
    xlen :
    ylen :
    hh :
    wd :
    space_srf :  min space between VM edge and SRF (km)
    space_land : min space between VM edge and land (km)
    logger :

    Returns
    -------
    origin, bearing, xlen, ylen : updated values

    """
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
            output=(wd / "tempEXT.tmp").as_posix(),
            min_edge_points=80,
            close=False,
        )
        isections, icomps = gmt.intersections(
            [
                (wd / "srf.path").as_posix(),
                NZ_LAND_OUTLINE.as_posix(),
                (wd / "tempEXT.tmp").as_posix(),
            ],
            items=True,
            containing=(wd / "tempEXT.tmp").as_posix(),
        )
        if len(isections) == 0:
            scan_extremes[i] = np.nan
            continue
        # shift different intersections by item-specific padding amount

        for c, intersection in enumerate(isections):
            if (wd / "srf.path").as_posix() in icomps[c]:
                diff = space_srf
            elif NZ_LAND_OUTLINE.as_posix() in icomps[c]:
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
            f"New origin after x shift: {origin}, new bearing after x shift: {bearing}"
        )

    xlen = math.ceil((dist_east_from_mid + dist_west_from_mid) / hh) * hh
    logger.debug(f"New xlen: {xlen}")

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
            f"New origin after y shift: {origin}, new bearing after y shift: {bearing}"
        )

    ylen -= over_s_max + over_n_max
    logger.debug(f"New ylen: {ylen}")

    return origin, bearing, xlen, ylen


def get_vm_land_proportion(
    c1: tuple, c2: tuple, c3: tuple, c4: tuple, wd: Path = Path(".")
):
    """
    Proportion of VM on land. Can be used to decide whether area optimisation should be done

    Parameters
    ----------
    c1 : tuple of (lon,lat)
    c2 : as above
    c3 : as above
    c4 : as above
    wd : working directory

    Returns
    -------
    Proportion (%)
    """
    path_vm = wd / "mask_vm.path"
    grd_vm = wd / "mask_vm.grd"
    grd_land = wd / "mask_land.grd"
    grd_vm_land = wd / "mask_vm_land.grd"

    geo.path_from_corners([c1, c2, c3, c4], output=path_vm)
    vm_region = corners2region(c1, c2, c3, c4)
    dxy = get_dxy(vm_region)

    gmt.grd_mask("f", grd_land, region=vm_region, dx=dxy, dy=dxy, wd=wd)
    gmt.grd_mask(path_vm, grd_vm, region=vm_region, dx=dxy, dy=dxy, wd=wd)
    gmt.grdmath([grd_land, grd_vm, "BITAND", "=", grd_vm_land], wd=wd)

    return grd_proportion(grd_vm_land) / grd_proportion(grd_vm) * 100


def centre_lon(lat_target: float):
    """
    get longitude on NZ centre path at given latitude

    Parameters
    ----------
    lat_target :

    Returns
    -------
    Interpolated longitude of a location whose latitiude is lat_target

    """

    lat = []
    lon = []

    with open(NZ_CENTRE_LINE, "r") as c:
        for line in c:
            ll = list(map(float, line.split()))
            lat.append(ll[1])
            lon.append(ll[0])

    lon_target = np.interp(lat_target, lat, lon)
    return lon_target


def optimise_vm_params(
    srf_meta: dict,
    ds_multiplier: float,
    dt: float,
    hh: float,
    pgv: float,
    temp_dir: Path,
    deep_rupture: bool = False,
    optimise: bool = True,
    target_land_coverage: float = 99.0,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Optimises VM domain and returns a dictionary that extends vm_params_dict

    Parameters
    ----------
    srf_meta : Data dict extracted from REL csv file
    ds_multiplier :
    dt :
    hh :
    pgv :
    deep_rupture : If true, continue even if the rupture is too deep.
    optimise : Performs area optimisation if set true
    target_land_coverage : Land coverage level (%) that triggers optimisation if not met. No guarantee to attain this level
    temp_dir : Temporary work directory (eg. .../Data/VMs/Hossack/_tmp_Hossack_nho8ui41)
    logger :

    Returns
    -------
    success: True/False to determine whether vm_params.yaml will be produced.
    vm_params_dict_extended: that contains extra info not stored in vm_params.yaml

    """

    # properties stored in classes (fault of external code)
    faultprop.Mw = srf_meta["mag"]
    faultprop.rake = srf_meta["rake"]
    faultprop.dip = srf_meta["dip"]
    faultprop.faultstyle = None

    success = True  # will produce VM

    # rrup to reach wanted PGV
    if pgv == -1.0:
        pgv = mag2pgv(faultprop.Mw)
        logger.debug(f"pgv determined from magnitude. pgv set to {pgv}")
    else:
        logger.debug(f"pgv set to {pgv}")
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
        rjb, hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=temp_dir
    )

    # for plotting and calculating VM domain distance
    with open(temp_dir / "srf.path", "wb") as sp:
        for plane in srf_meta["corners"]:
            sp.write("> srf plane\n".encode())
            np.savetxt(sp, plane, fmt="%f")
            sp.write("%f %f\n".encode() % (tuple(plane[0])))

    if fault_depth > rrup and not deep_rupture:
        logger.warning(
            "fault_depth > rrup. Fault too deep, will not generate VM. "
            f"Setting xlen, ylen to 0. Previous values: x:{xlen0}, y:{ylen0}"
        )
        xlen0 = 0
        ylen0 = 0
        success = False

    o1, o2, o3, o4 = geo.build_corners(origin, bearing, xlen0, ylen0)
    vm0_region = corners2region(o1, o2, o3, o4)
    plot_region = (
        vm0_region[0] - 1,
        vm0_region[1] + 1,
        vm0_region[2] - 1,
        vm0_region[3] + 1,
    )

    # proportion in ocean
    if xlen0 <= 0 or ylen0 <= 0:
        logger.debug(f"Optimising skipped for {srf_meta['name']}. 100% Land coverage")
        land0 = 100
        optimise = False
    else:
        land0 = get_vm_land_proportion(o1, o2, o3, o4, wd=temp_dir)
        logger.debug(f"Land coverage found to be {land0}%")

    if not optimise:
        logger.info(f"No optimisation for {srf_meta['name']} : no_optimise==True")

    if faultprop.Mw < 3.5:
        optimise = False
        logger.info(f"No optimisation for {srf_meta['name']} : mw<3.5")

    if land0 >= target_land_coverage:
        logger.info(
            f"No optimisation for {srf_meta['name']} : land coverage >= {target_land_coverage}%"
        )
        optimise = False

    # modify VM if necessary
    if optimise:
        logger.info(
            f"Optimising : {srf_meta['name']} has mw>= 3.5 and land coverage < {target_land_coverage}%"
        )

        # rotation based on centre line bearing at this point
        l1 = centre_lon(vm0_region[2])
        l2 = centre_lon(vm0_region[3])
        mid = geo.ll_mid(l1, vm0_region[2], l2, vm0_region[3])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, vm0_region[3]))

        # wanted distance is at corners, not middle top to bottom
        _, xlen1, ylen1 = determine_vm_extent(
            rjb, hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=temp_dir
        )

        logger.debug(
            f"Pre-optimisation origin: {origin}, bearing: {bearing}, xlen: {xlen1}, ylen: {ylen1}"
        )
        # cut down ocean areas
        (origin, bearing, xlen1, ylen1) = reduce_domain(
            origin, bearing, xlen1, ylen1, hh, temp_dir, logger=logger
        )
        logger.debug(
            f"After optimisation origin: {origin}, bearing: {bearing}, xlen: {xlen1}, ylen: {ylen1}"
        )

        try:
            c1, c2, c3, c4 = geo.build_corners(origin, bearing, xlen1, ylen1)
        except ValueError:
            error_message = f"Error for vm {srf_meta['name']}. Raising exception."
            logger.log(qclogging.NOPRINTCRITICAL, error_message)
            raise ValueError(error_message)
        else:
            logger.debug(f"Optimised corners of the VM are: {(c1, c2, c3, c4)}")

        # proportion in ocean
        land1 = get_vm_land_proportion(c1, c2, c3, c4, wd=temp_dir)
        logger.debug(f"Optimised land coverage found to be {land1}%")

        # adjust region to fit new corners
        plot_region = (
            min(c4[0], c1[0], c3[0], c2[0], plot_region[0]),
            max(c4[0], c1[0], c3[0], c2[0], plot_region[1]),
            min(c4[1], c1[1], c3[1], c2[1], plot_region[2]),
            max(c4[1], c1[1], c3[1], c2[1], plot_region[3]),
        )
    else:  # not optimised
        # store original variables as final versions
        ylen1 = ylen0
        xlen1 = xlen0
        land1 = land0
        c1, c2, c3, c4 = o1, o2, o3, o4

    # not enough land in final domain
    if math.floor(land1) == 0:
        logger.info(
            "Land coverage is less than 1%. Setting xlen and ylen to 0. Not creating VM"
        )
        xlen1 = 0
        ylen1 = 0
        success = False

    # zlen is independent from xlen and ylen
    zlen = round(get_max_depth(faultprop.Mw, srf_meta["dbottom"]) / hh) * hh
    logger.debug(f"zlen set to {zlen}")

    if xlen1 == 0 or ylen1 == 0 or zlen == 0:
        logger.debug(
            f"All xlen={xlen1} ylen={ylen1} zlen={zlen} should be non-zero. Not creating VM"
        )
        success = False

    bounds_valid = validate_vm_bounds([c1, c2, c3, c4], srf_meta["corners"].tolist())
    if not bounds_valid:
        logger.warning(f"Bounds not valid, not making VM: {bounds_valid}")
        success = False

    # modified sim time
    vm_corners = np.asarray([c1, c2, c3, c4])
    initial_time = get_sim_duration(
        vm_corners,
        np.concatenate(srf_meta["corners"], axis=0),
        ds_multiplier,
        fault_depth,
        logger=logger,
    )
    sim_duration = (initial_time // dt) * dt
    logger.debug(
        f"Unrounded sim duration: {initial_time}. Rounded sim duration: {sim_duration}"
    )

    # optimisation results
    return (
        success,
        {
            "name": srf_meta["name"],
            "mag": faultprop.Mw,
            "dbottom": srf_meta["dbottom"],
            "zlen": zlen,
            "sim_time": sim_duration,
            "xlen": xlen0,
            "ylen": ylen0,
            "land": land0,
            "zlen_mod": zlen,
            "sim_time_mod": sim_duration,
            "xlen_mod": xlen1,
            "ylen_mod": ylen1,
            "land_mod": land1,
            "origin": origin,
            "bearing": bearing,
            "adjusted": optimise,
            "plot_region": plot_region,
            "path": "{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n".format(
                o1[0], o1[1], o2[0], o2[1], o3[0], o3[1], o4[0], o4[1]
            ),
            "path_mod": "{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n".format(
                c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]
            ),
        },
    )


def main(
    srf_meta: dict,
    ds_multiplier: float,
    dt: float,
    hh: float,
    min_vs: float,
    outdir: Path,
    pgv: float,
    vm_topo: str,
    vm_version: str,
    deep_rupture: bool = False,
    target_land_coverage: float = 99.0,
    optimise: bool = True,
    plot_enabled: bool = True,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Orchestrates conversion from rel CSV file to vm_params.yaml

    Parameters
    ----------
    srf_meta : Data extracted from REL csv file
    ds_multiplier :
    dt :
    hh :
    min_vs :
    outdir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    pgv :
    vm_topo :
    vm_version :
    deep_rupture :  If true, continue even if the rupture is too deep.
    target_land_coverage : Land coverage level (%) that triggers optimisation if not met.
    optimise : Performs area optimisation if set true
    plot_enabled : If true, plot the vm domain
    logger:
    """

    # temp directory for current process
    with TemporaryDirectory(prefix=f"_tmp_{srf_meta['name']}_", dir=outdir) as temp_dir:
        temp_dir = Path(temp_dir)

        qclogging.add_general_file_handler(
            logger, outdir / f"rel2vm_params_{srf_meta['name']}_log.txt"
        )
        vm_params_path = outdir / "vm_params.yaml"

        success, vm_params_dict_extended = optimise_vm_params(
            srf_meta,
            ds_multiplier,
            dt,
            hh,
            pgv,
            temp_dir,
            deep_rupture=deep_rupture,
            optimise=optimise,
            target_land_coverage=target_land_coverage,
            logger=logger,
        )

        if success:
            xlen = float(vm_params_dict_extended["xlen_mod"])
            ylen = float(vm_params_dict_extended["ylen_mod"])
            zmax = float(vm_params_dict_extended["zlen_mod"])
            zmin = 0.0

            code = "rt"
            vm_params_dict = {
                "mag": float(faultprop.Mw),
                "MODEL_LAT": float(vm_params_dict_extended["origin"][1]),
                "MODEL_LON": float(vm_params_dict_extended["origin"][0]),
                "MODEL_ROT": float(vm_params_dict_extended["bearing"]),
                "hh": hh,
                "min_vs": min_vs,
                "model_version": vm_version,
                "topo_type": vm_topo,
                "extent_x": xlen,
                "extent_y": ylen,
                "extent_zmax": zmax,
                "extent_zmin": zmin,
                "sim_duration": float(vm_params_dict_extended["sim_time_mod"]),
                "flo": min_vs / (5.0 * hh),
                "nx": int(round(float(xlen) / hh)),
                "ny": int(round(float(ylen) / hh)),
                "nz": int(round(float(zmax - zmin) / hh)),
                "sufx": f"_{code}01-h{hh:.3f}",
            }

            if vm_params_path.exists():
                logger.info(
                    f"Warning: {vm_params_path} already exists and is overwritten."
                )
                os.remove(vm_params_path)

            dump_yaml(vm_params_dict, vm_params_path)
            logger.debug(f"Saved {vm_params_path}")
            print(f"Success: Wrote vm_params.yaml at {outdir}", file=sys.stderr)

            # generate a corners like NZVM would have
            logger.debug("Saving VeloModCorners.txt")
            write_corner_file(outdir, vm_params_dict_extended["path_mod"])
        else:
            logger.error("Failed: Not good VM params to proceed")

        if plot_enabled:
            plot_vm(
                vm_params_dict_extended,
                srf_meta["corners"],
                NZ_LAND_OUTLINE,
                NZ_CENTRE_LINE,
                faultprop.Mw,
                outdir,
                temp_dir,
                logger=logger,
            )


def write_corner_file(outdir: Path, paths: str):
    """
    Write a corner file
    Parameters
    ----------
    outdir :
    paths : str of corner coordinates "Lon1\tLat1\nLon2\tLat2\nLon3\tLat3\nLon4\tLat4" where Lat or Lon are in .6f format.
    """
    with open(f"{outdir}/VeloModCorners.txt", "w") as c:
        c.write(">Velocity model corners(python generated)\n")
        c.write(">Lon\tLat\n")
        c.write(paths)


def load_rel(rel_file: Path, logger: Logger = qclogging.get_basic_logger()):
    """

    Parameters
    ----------
    rel_file : Path to REL csv file
    logger :

    Returns
    -------
    SRF meta data dictionary
    """

    # name of the fault is found from the basename

    name = get_fault_from_realisation(rel_file)  # XXXX_REL01.csv --> XXXX
    logger.debug(f"Found first REL file for {name}")

    rel_df = pd.read_csv(rel_file)

    # common attributes in all types of rel csv
    def rel_meta(attr):
        return rel_df[attr].loc[0] if attr in rel_df.columns else None

    type = rel_meta("type")

    if type == 1:  # point source
        hdepth = rel_meta("depth")
        planes = [
            {
                "centre": [rel_meta("longitude"), rel_meta("latitude")],
                "length": 0.1,
                "width": 0.1,
                "strike": rel_meta("strike"),
                "dip": rel_meta("dip"),
                "dtop": hdepth,
            }
        ]
        corners, [dbottom] = get_corners_dbottom(planes)

    elif type == 2:  # finite fault single plane
        planes = [
            {
                "centre": [rel_meta("longitude"), rel_meta("latitude")],
                "length": rel_meta("flen"),
                "width": rel_meta("fwid"),
                "strike": rel_meta("strike"),
                "dip": rel_meta("dip"),
                "dtop": rel_meta("dtop"),
            }
        ]
        corners, [dbottom] = get_corners_dbottom(planes)

    else:  # type 4 (we never have type 3)
        assert type != 3

        planes = []
        num_subfaults = rel_meta("plane_count")
        for i in range(num_subfaults):
            plane = {}
            plane["centre"] = [
                rel_meta(f"clon_subfault_{i}"),
                rel_meta(f"clat_subfault_{i}"),
            ]
            plane["length"] = rel_meta(f"length_subfault_{i}")
            plane["width"] = rel_meta(f"width_subfault_{i}")
            plane["strike"] = rel_meta(f"strike_subfault_{i}")
            plane["dip"] = rel_meta(f"dip_subfault_{i}")
            plane["dtop"] = rel_meta(f"dtop_subfault_{i}")
            planes.append(plane)

        corners, _ = get_corners_dbottom(planes, dip_dir=rel_meta("dip_dir"))
        dbottom = rel_meta("dbottom")

    if type >= 2:  # type 2 and 4 will have hdepth computed
        hdepth = np.round(rel_meta("dhypo"), decimals=1) * np.sin(
            np.radians(rel_meta("dip"))
        ) + rel_meta("dtop")

    return {
        "name": name,
        "dip": rel_meta("dip"),
        "rake": rel_meta("rake"),
        "dbottom": dbottom,
        "corners": np.array(corners),
        "mag": rel_meta("magnitude"),
        "hdepth": hdepth,
    }


def load_args():
    """
    Unpacks arguments and does basic checks

    Parameters
    ----------

    Returns
    -------
    Processed arguments

    """
    parser = ArgumentParser()
    arg = parser.add_argument

    arg("rel_file", help="REL csv file")

    arg(
        "-o",
        "--outdir",
        help="output directory to place VM files "
        "(if not specified, the same location as rel_file is in)",
        default=None,
    )

    arg(
        "--pgv",
        help="max PGV at velocity model perimiter (estimated, cm/s)",
        type=float,
        default=-1.0,
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
        choices=[x.value for x in constants.VelocityModelTopographyType],
        default="BULLDOZED",
    )
    arg(
        "--no-optimise",
        help="Don't try and optimise the vm if it is off shore. Removes dependency on having GMT coastline data",
        action="store_false",
        default=True,
        dest="optimise",
    )
    arg(
        "--deep-rupture",
        help="Continue even if too deep",
        action="store_true",
        default=False,
    )
    arg(
        "--target-land-coverage",
        help="Land coverage level (%) that triggers optimisation if not met (Default: 99.0)",
        type=float,
        default=None,
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

    args.rel_file = Path(args.rel_file).resolve()
    if not args.rel_file.exists():
        parser.error(f"REL csv file not found: {args.rel_file}")

    if args.outdir is None:
        args.outdir = args.rel_file.parent
    args.outdir = Path(args.outdir).resolve()

    if args.optimise is False and args.target_land_coverage is not None:
        parser.error(
            "Both no_optimise and target_land_coverage are set. Set either one"
        )

    if args.target_land_coverage is None:
        args.target_land_coverage = 99.0  # default if not set

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("rel2vm_params")
    qclogging.add_general_file_handler(logger, Path.cwd() / "rel2vm_params_log.txt")

    args = load_args()

    logger.debug("Loading REL csv")
    srf_meta = load_rel(args.rel_file, logger=logger)

    # prepare to run
    os.makedirs(args.outdir, exist_ok=True)

    main(
        srf_meta,
        args.ds_multiplier,
        args.dt,
        args.hh,
        args.min_vs,
        args.outdir,
        args.pgv,
        args.vm_topo,
        args.vm_version,
        deep_rupture=args.deep_rupture,
        optimise=args.optimise,
        logger=logger,
    )

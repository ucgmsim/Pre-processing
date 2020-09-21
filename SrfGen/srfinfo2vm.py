#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example (set out_dir with -o or --out-dir):
mpirun -n 3 srfinfo2vm.py "Srf/*.info"

HELP:
./srfinfo2vm -h
"""

from argparse import ArgumentParser
from distutils.spawn import find_executable
from glob import glob
import math
from logging import Logger
from multiprocessing import Pool
import os
import platform
import subprocess
from shutil import rmtree, move
from subprocess import Popen
import sys
from tempfile import mkdtemp

from h5py import File as h5open
import numpy as np

from qcore import constants, geo, gmt, qclogging
from qcore.geo import R_EARTH
from qcore.utils import dump_yaml
from qcore.validate_vm import validate_vm

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

from createSRF import leonard, mag2mom, mom2mag
from gen_coords import gen_coords

script_dir = os.path.dirname(os.path.abspath(__file__))
NZ_CENTRE_LINE = os.path.join(script_dir, "NHM/res/centre.txt")
NZ_LAND_OUTLINE = os.path.join(script_dir, "NHM/res/rough_land.txt")
NZVM_BIN = find_executable("NZVM")

siteprop = Site()
siteprop.vs30 = 500

faultprop = Fault()
# DON'T CHANGE THIS - this assumes we have a surface point source
faultprop.ztor = 0.0

S_WAVE_KM_PER_S = 3.5


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


def save_vm_config(
    nzvm_cfg=None,
    vm_params=None,
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
    Store VM config for NZVM generator and vm_params metadata store.
    nzvm_cfg: path to NZVM cfg file
    vm_params: path to vm_params.py that stores metadata
    vm_dir: folder for NZVM output (cannot exist)
    origin: model origin (longitude, latitude)
    rot: model rotation
    """
    assert vm_dir is not None
    if nzvm_cfg is not None:
        assert not os.path.exists(vm_dir)
        with open(nzvm_cfg, "w") as vmd:
            vmd.write(
                "\n".join(
                    [
                        "CALL_TYPE=GENERATE_VELOCITY_MOD",
                        "MODEL_VERSION=%s" % (model_version),
                        "OUTPUT_DIR=%s" % (vm_dir),
                        "ORIGIN_LAT=%s" % (origin[1]),
                        "ORIGIN_LON=%s" % (origin[0]),
                        "ORIGIN_ROT=%s" % (rot),
                        "EXTENT_X=%s" % (xlen),
                        "EXTENT_Y=%s" % (ylen),
                        "EXTENT_ZMAX=%s" % (zmax),
                        "EXTENT_ZMIN=%s" % (zmin),
                        "EXTENT_Z_SPACING=%s" % (hh),
                        "EXTENT_LATLON_SPACING=%s" % (hh),
                        "MIN_VS=%s" % (min_vs),
                        "TOPO_TYPE=%s\n" % (topo_type),
                    ]
                )
            )
        logger.debug("Saved nzvm config file to {}".format(nzvm_cfg))
    if vm_params is not None:
        # must also convert mag from np.float to float
        dump_yaml(
            {
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
            },
            "{}.yaml".format(vm_params),
        )
        logger.debug("Saved vm_params.yaml to {}".format("{}.yaml".format(vm_params)))


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


# maximum width along lines a and b
def reduce_domain(
    origin,
    bearing,
    xlen,
    ylen,
    hh,
    space_srf,
    space_land,
    wd,
    logger: Logger = qclogging.get_basic_logger(),
):
    """Reduces the domain of the VM by removing areas that are over sea only. Gives a buffer around coast line and the
    srf boundary. Returns the new values to be used.
    Uses great circle geometry """
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


def gen_vm(
    args,
    srf_meta,
    vm_params_dict,
    mag,
    ptemp,
    logger: Logger = qclogging.get_basic_logger(),
):
    # store configs
    out_vm_dir = os.path.join(args.out_dir, srf_meta["name"])
    logger.info("Generating VM. Saving it to {}".format(out_vm_dir))
    vm_params_dict["vm_dir"] = out_vm_dir
    vm_working_dir = os.path.join(out_vm_dir, "output")
    os.makedirs(vm_working_dir, exist_ok=True)
    nzvm_cfg = os.path.join(ptemp, "nzvm.cfg")
    vm_params_path = os.path.join(ptemp, "vm_params")
    # NZVM won't run if folder exists
    if os.path.exists(vm_working_dir):
        logger.debug("VM working directory {} already exists.".format(vm_working_dir))
        rmtree(vm_working_dir)
    save_vm_config(
        nzvm_cfg=nzvm_cfg,
        vm_params=vm_params_path,
        vm_dir=vm_working_dir,
        origin=vm_params_dict["origin"],
        rot=vm_params_dict["bearing"],
        xlen=vm_params_dict["xlen_mod"],
        ylen=vm_params_dict["ylen_mod"],
        zmax=vm_params_dict["zlen_mod"],
        hh=args.hh,
        min_vs=args.min_vs,
        mag=mag,
        centroid_depth=srf_meta["hdepth"],
        sim_duration=vm_params_dict["sim_time_mod"],
        topo_type=args.vm_topo,
        model_version=args.vm_version,
    )
    if args.novm:
        logger.debug(
            "--novm set, generating configuration files, but not generating VM."
        )
        # save important files
        logger.debug("Creating directory {}".format(vm_working_dir))
        os.makedirs(vm_working_dir)
        move(nzvm_cfg, vm_working_dir)
        move("{}.yaml".format(vm_params_path), vm_working_dir)
        logger.debug(
            "Moved nvzm config and vm_params yaml to {}".format(vm_working_dir)
        )
        # generate a corners like NZVM would have
        logger.debug("Saving VeloModCorners.txt")
        with open("{}/VeloModCorners.txt".format(vm_params_dict["vm_dir"]), "wb") as c:
            c.write("> VM corners (python generated)\n".encode())
            c.write(">Lon    Lat\n".encode())
            c.write(vm_params_dict["path_mod"].encode())
        return

    # NZVM won't find resources if WD is not NZVM dir, stdout not MPROC friendly
    with open(os.path.join(ptemp, "NZVM.out"), "w") as logfile:
        logger.debug("Running NZVM binary")
        nzvm_env = os.environ.copy()
        nzvm_env["OMP_NUM_THREADS"] = str(args.vm_threads)
        nzvm_exe = Popen(
            [NZVM_BIN, nzvm_cfg],
            cwd=os.path.dirname(NZVM_BIN),
            stdout=logfile,
            env=nzvm_env,
        )
        nzvm_exe.communicate()
    logger.debug("Moving VM files to vm directory")
    # fix up directory contents
    move(os.path.join(vm_working_dir, "Velocity_Model", "rho3dfile.d"), out_vm_dir)
    move(os.path.join(vm_working_dir, "Velocity_Model", "vp3dfile.p"), out_vm_dir)
    move(os.path.join(vm_working_dir, "Velocity_Model", "vs3dfile.s"), out_vm_dir)
    move(os.path.join(vm_working_dir, "Log", "VeloModCorners.txt"), out_vm_dir)
    logger.debug("Removing Log and Velocity_Model directories")
    logger.debug("Moving nzvm config and vm_params yaml to vm directory")
    move(nzvm_cfg, out_vm_dir)
    move("%s.yaml" % (vm_params_path), out_vm_dir)
    rmtree(vm_working_dir)
    # create model_coords, model_bounds etc...
    logger.debug("Generating coords")
    gen_coords(vm_dir=out_vm_dir)
    # validate
    logger.debug("Validating vm")
    success, message = validate_vm(out_vm_dir)
    if success:
        vm_check_str = f"VM check passed: {vm_working_dir}"
        logger.debug(vm_check_str)
        sys.stderr.write(vm_check_str)
    else:
        logger.log(
            qclogging.NOPRINTCRITICAL,
            "VM check for {} failed: {}".format(vm_working_dir, message),
        )
        sys.stderr.write("VM check BAD: {}\n".format(message))


def plot_vm(
    vm_params, srf_corners, mag, ptemp, logger: Logger = qclogging.get_basic_logger()
):
    logger.debug("Plotting vm")
    p = gmt.GMTPlot(os.path.join(ptemp, "optimisation.ps"))
    p.spacial("M", vm_params["plot_region"], sizing=7)
    p.coastlines()

    # filled slip area
    p.path("%s/srf.path" % (ptemp), is_file=True, fill="yellow", split="-")
    # top edge
    for plane in srf_corners:
        p.path("\n".join([" ".join(map(str, ll)) for ll in plane[:2]]), is_file=False)

    # vm domain (simple and adjusted)
    p.path(vm_params["path"], is_file=False, close=True, fill="black@95")
    if vm_params["adjusted"]:
        p.path(
            vm_params["path_mod"],
            is_file=False,
            close=True,
            fill="black@95",
            split="-",
            width="1.0p",
        )

    # info text for simple and adjusted domains
    p.text(
        sum(vm_params["plot_region"][0:2]) / 2.0,
        vm_params["plot_region"][3],
        "Mw: %.2f X: %.0fkm, Y: %.0fkm, land: %.0f%%"
        % (mag, vm_params["xlen"], vm_params["ylen"], vm_params["land"]),
        align="CT",
        dy=-0.1,
        box_fill="white@50",
    )
    if vm_params["adjusted"]:
        p.text(
            sum(vm_params["plot_region"][0:2]) / 2.0,
            vm_params["plot_region"][3],
            "MODIFIED land: %.0f%%" % (vm_params["land_mod"]),
            align="CT",
            dy=-0.25,
            box_fill="white@50",
        )

    # land outlines blue, nz centre line (for bearing calculation) red
    p.path(NZ_LAND_OUTLINE, is_file=True, close=False, colour="blue", width="0.2p")
    p.path(NZ_CENTRE_LINE, is_file=True, close=False, colour="red", width="0.2p")

    # actual corners retrieved from NZVM output or generated if args.novm
    # not available if VM was skipped
    vm_exists = "vm_dir" in vm_params
    if vm_exists:
        logger.debug("vm exists, getting corners from VeloModCorners")
        p.points(
            "%s/VeloModCorners.txt" % (vm_params["vm_dir"]),
            fill="red",
            line=None,
            shape="c",
            size=0.05,
        )
    else:
        logger.debug("vm doesn't exist, deriving corners from path mod")
        p.points(
            vm_params["path_mod"],
            is_file=False,
            fill="red",
            line=None,
            shape="c",
            size="0.05",
        )

    # store PNG
    p.finalise()
    logger.debug("Saving image")
    if vm_exists:
        p.png(dpi=200, clip=True, background="white", out_dir=vm_params["vm_dir"])
    else:
        p.png(
            dpi=200,
            clip=True,
            background="white",
            out_name=os.path.abspath(os.path.join(ptemp, os.pardir, vm_params["name"])),
        )


def create_vm(args, srf_meta, logger_name: str = "srfinfo2vm"):
    # temp directory for current process
    logger = qclogging.get_realisation_logger(
        qclogging.get_logger(logger_name), srf_meta["name"]
    )
    ptemp = mkdtemp(prefix="_tmp_%s_" % (srf_meta["name"]), dir=args.out_dir)

    # properties stored in classes (fault of external code)
    faultprop.Mw = srf_meta["mag"]
    faultprop.rake = srf_meta["rake"]
    faultprop.dip = srf_meta["dip"]
    faultprop.faultstyle = None
    # rrup to reach wanted PGV
    if args.pgv == -1.0:
        pgv = mag2pgv(faultprop.Mw)
        logger.debug("pgv determined from magnitude. pgv set to {}".format(pgv))
    else:
        pgv = args.pgv
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
            args.min_rjb, rrup
        )  # sets rrup equal to rjb to ensure deep ruptures have sufficient VM size

    # original, unrotated vm
    bearing = 0
    origin, xlen0, ylen0 = determine_vm_extent(
        rjb, args.hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=ptemp
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
    if args.no_optimise or xlen0 <= 0 or ylen0 <= 0:
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
            rjb, args.hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=ptemp
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
            args.hh,
            args.space_srf,
            args.space_land,
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
    zlen = round(auto_z(faultprop.Mw, srf_meta["dbottom"]) / args.hh) * args.hh
    logger.debug("zlen set to {}".format(zlen))
    # modified sim time
    vm_corners = np.asarray([c1, c2, c3, c4])
    initial_time = auto_time2(
        vm_corners,
        np.concatenate(srf_meta["corners"], axis=0),
        args.ds_multiplier,
        fault_depth,
        logger=logger,
    )
    sim_time1 = (initial_time // args.dt) * args.dt
    logger.debug(
        "Unrounded sim duration: {}. Rounded sim duration: {}".format(
            initial_time, sim_time1
        )
    )

    # optimisation results
    vm_params = {
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

    # will not create VM if it is entirely in the ocean
    if xlen1 != 0 and ylen1 != 0 and zlen != 0:
        # run the actual generation
        gen_vm(args, srf_meta, vm_params, faultprop.Mw, ptemp, logger=logger)
    else:
        logger.debug("At least one dimension was 0. Not generating VM")
    # plot results
    plot_vm(vm_params, srf_meta["corners"], faultprop.Mw, ptemp, logger=logger)

    # working dir cleanup, return info about VM
    logger.debug("Cleaning up temp directory {}".format(ptemp))
    rmtree(ptemp)
    return vm_params


def load_msgs(args, logger: Logger = qclogging.get_basic_logger()):
    # returns list of appropriate srf metadata
    msgs = []
    faults = set()

    # add task for every info file
    info_files = glob(args.info_glob)
    for info in info_files:
        if not os.path.exists(info):
            logger.info("SRF info file not found: {}".format(info))
            continue

        # name is unique and based on basename
        name = os.path.splitext(os.path.basename(info))[0].split("_")[0]

        if name in faults:
            continue
        logger.debug("Found first info file for {}".format(name))
        faults.add(name)

        with h5open(info, "r") as h:
            a = h.attrs
            try:
                rake = a["rake"][0][0]
            except (TypeError, IndexError):
                rake = a["rake"]
            try:
                mag = mom2mag(sum(map(mag2mom, a["mag"])))
            except TypeError:
                mag = a["mag"]
            msgs.append(
                (
                    args,
                    {
                        "name": name,
                        "dip": a["dip"][0],
                        "rake": rake,
                        "dbottom": a["dbottom"][0],
                        "corners": a["corners"],
                        "mag": mag,
                        "hdepth": a["hdepth"],
                    },
                    qclogging.get_realisation_logger(logger, name).name,
                )
            )
    return msgs


def load_msgs_nhm(args, logger: Logger = qclogging.get_basic_logger()):
    msgs = []
    logger.debug("Loading NHM file {}".format(args.nhm_file))
    with open(args.nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
    n_i = 0
    while n_i < len(nhm):
        n = nhm[n_i].strip()
        logger.debug("Getting data for {}".format(n))
        dip = float(nhm[n_i + 3].split()[0])
        dip_dir = float(nhm[n_i + 4])
        rake = float(nhm[n_i + 5])
        dbottom = float(nhm[n_i + 6].split()[0])
        dtop = float(nhm[n_i + 7].split()[0])
        n_pt = int(nhm[n_i + 11])
        trace = nhm[n_i + 12 : n_i + 12 + n_pt]
        trace = [list(map(float, pair)) for pair in map(str.split, trace)]

        # derived properties
        dbottom += 3 * (dbottom >= 12)
        pwid = (dbottom - dtop) / math.tan(math.radians(dip))
        fwid = (dbottom - dtop) / math.sin(math.radians(dip))
        trace_length = sum(
            [
                geo.ll_dist(trace[i][0], trace[i][1], trace[i + 1][0], trace[i + 1][1])
                for i in range(n_pt - 1)
            ]
        )
        corners = np.zeros((n_pt - 1, 4, 2))
        for i in range(len(corners)):
            corners[i, :2] = trace[i : i + 2]
            corners[i, 2] = geo.ll_shift(
                corners[i, 1, 1], corners[i, 1, 0], pwid, dip_dir
            )[::-1]
            corners[i, 3] = geo.ll_shift(
                corners[i, 0, 1], corners[i, 0, 0], pwid, dip_dir
            )[::-1]

        # move to next fault
        n_i += 13 + n_pt

        mag = leonard(rake, fwid * trace_length)
        msgs.append(
            (
                args,
                {
                    "name": n,
                    "dip": dip,
                    "rake": rake,
                    "dbottom": dbottom,
                    "corners": corners,
                    "mag": mag,
                    "hdepth": dtop + 0.5 * (dbottom - dtop),
                },
                qclogging.get_realisation_logger(logger, n).name,
            )
        )

    return msgs


def store_nhm_selection(
    out_dir, reports, logger: Logger = qclogging.get_basic_logger()
):
    selected = os.path.join(out_dir, "nhm_selection.txt")
    excluded = os.path.join(out_dir, "excluded.txt")
    logger.debug("Saving selection and exclusion files")
    with open(selected, "w") as sf:
        with open(excluded, "w") as ef:
            for r in reports:
                if r["xlen_mod"] == 0 or r["ylen_mod"] == 0 or r["zlen"] == 0:
                    logger.debug(
                        "{} has atleast one dimension of length 0, excluding".format(
                            r["name"]
                        )
                    )
                    ef.write("%s\n" % (r["name"]))
                else:
                    logger.debug(
                        "{} has length greater than 0 for all dimensions, including".format(
                            r["name"]
                        )
                    )
                    sf.write(
                        "%s %dr\n" % (r["name"], min(max(r["mag"] * 20 - 110, 10), 50))
                    )
    logger.debug("Saving nhm selection file to {}".format(selected))
    logger.debug("Saving nhm exclusion file to {}".format(excluded))


def store_summary(table, info_store, logger: Logger = qclogging.get_basic_logger()):
    # initialise table file
    logger.debug("Saving summary")
    with open(table, "w") as t:
        t.write(
            '"name","mw","plane depth (km, to bottom)",'
            '"vm depth (zlen, km)","xlen (km)",'
            '"ylen (km)","land (% cover)","adjusted vm depth (zlen, km)",'
            '"adjusted sim time (s)","adjusted xlen (km)",'
            '"adjusted ylen (km)","adjusted land (% cover)"\n'
        )

        for i in info_store:
            t.write(
                "%s,%s,%s,%s,%s,%s,%.0f,%s,%s,%s,%s,%.0f\n"
                % (
                    i["name"],
                    i["mag"],
                    i["dbottom"],
                    i["zlen"],
                    i["xlen"],
                    i["ylen"],
                    i["land"],
                    i["zlen_mod"],
                    i["sim_time_mod"],
                    i["xlen_mod"],
                    i["ylen_mod"],
                    i["land_mod"],
                )
            )
    logger.info("Saved summary to {}".format(table))


def load_args(parser=ArgumentParser(), logger: Logger = qclogging.get_basic_logger()):

    arg = parser.add_argument
    arg("info_glob", help="info file selection expression. eg: Srf/*.info")
    arg(
        "--nhm-file",
        help="path to NHM if using info_glob == 'NHM'",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "NHM", "NZ_FLTmodel_2010.txt"
        ),
    )
    parser.add_argument(
        "-o", "--out-dir", help="directory to place outputs", default="VMs"
    )
    arg(
        "--pgv",
        help="max PGV at velocity model perimiter (estimated, cm/s)",
        type=float,
        default=5.0,
    )
    arg("--hh", help="velocity model grid spacing (km)", type=float, default=0.4)
    arg(
        "--dt",
        help="timestep to estimate simulation duration (s)",
        type=float,
        default=0.005,
    )
    arg(
        "--space-land",
        help="min space between VM edge and land (km)",
        type=float,
        default=5.0,
    )
    arg(
        "--space-srf",
        help="min space between VM edge and SRF (km)",
        type=float,
        default=15.0,
    )
    arg("--min-vs", help="for nzvm gen and flo (km/s)", type=float, default=0.5)
    arg("-n", "--nproc", help="number of processes", type=int, default=1)
    arg(
        "-t",
        "--vm_threads",
        "--threads",
        help="number of threads for the VM generation",
        type=int,
        default=1,
    )
    arg("--novm", help="only generate parameters", action="store_true")
    arg("--vm-version", help="velocity model version to generate", default="1.65")
    arg(
        "--vm-topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )
    arg("--selection", help="also generate NHM selection file", action="store_true")
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
    args.out_dir = os.path.abspath(args.out_dir)

    if not args.novm:
        if NZVM_BIN is None:
            message = """NZVM binary not in PATH
You can compile it from here: https://github.com/ucgmsim/Velocity-Model
Then add it to the PATH by either:
sudo cp /location/to/Velocity-Model/NZVM /usr/bin/
or:
export PATH=$PATH:/location/to/Velocity-Model"""
            logger.log(qclogging.NOPRINTERROR, message)
            parser.error(message)

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("srfinfo2vm")
    qclogging.add_general_file_handler(
        logger, os.path.join(os.getcwd(), "srfinfo2vm_log.txt")
    )
    args = load_args(logger=logger)

    # load wanted fault information
    if args.info_glob == "NHM":
        logger.debug(
            "info_glob is NHM. Loading messages from NHM file: {}".format(args.nhm_file)
        )
        msg_list = load_msgs_nhm(args, logger=logger)
    else:
        logger.debug(
            "info_glob is not NHM, assuming it is a valid path (possibly with stars)"
        )
        msg_list = load_msgs(args, logger=logger)
    if len(msg_list) == 0:
        message = "Found nothing to do, exiting."
        logger.log(qclogging.NOPRINTCRITICAL, message)
        sys.exit(message)

    # prepare to run
    if not os.path.isdir(args.out_dir):
        logger.debug(
            "Output directory {} does not exist. Creating it now.".format(args.out_dir)
        )
        os.makedirs(args.out_dir)

    # distribute work
    logger.debug(
        "Number of processes to use given as {}, number of tasks: {}.".format(
            args.nproc, len(msg_list)
        )
    )
    p = Pool(processes=args.nproc)
    reports = p.starmap(create_vm, msg_list)

    # nhm selection formatted file and list of excluded VMs
    if args.selection:
        logger.debug("Storing NHM selection file")
        store_nhm_selection(args.out_dir, reports, logger=logger)
    # store summary
    store_summary(os.path.join(args.out_dir, "vminfo.csv"), reports, logger=logger)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)

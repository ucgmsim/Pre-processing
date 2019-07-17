#!/usr/bin/env python
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
from multiprocessing import Pool
import os
import platform
import subprocess
from shutil import rmtree, move
from subprocess import Popen, PIPE
import sys
from tempfile import mkdtemp
from time import time

from h5py import File as h5open
import numpy as np

from createSRF import leonard, mag2mom, mom2mag

# qcore library should already be in path
try:
    from qcore import constants
    from qcore import geo
    from qcore import gmt
    from gen_coords import gen_coords
    from qcore.validate_vm import validate_vm
    from qcore.utils import dump_yaml
except ImportError:
    sys.exit("""
qcore library has not been installed or is an old version.
you can install it using setup.py or add qcore/qcore to the PYTHONPATH like:
$ export PYTHONPATH=$PYTHONPATH:/location/to/qcore/qcore
qcore is available at https://github.com/ucgmsim/qcore""")

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

script_dir = os.path.dirname(os.path.abspath(__file__))
NZ_CENTRE_LINE = os.path.join(script_dir, "NHM/res/centre.txt")
NZ_LAND_OUTLINE = os.path.join(script_dir, "NHM/res/rough_land.txt")
NZVM_BIN = find_executable("NZVM")

siteprop = Site()
siteprop.vs30 = 500

faultprop = Fault()
# TODO ztor should be read from srfinfo file
faultprop.ztor = 0.0

# default scaling relationship
def mag2pgv(mag):
    return np.interp(
        mag, [3.5, 4.1, 4.7, 5.2, 6.4, 7.5], [0.015, 0.0375, 0.075, 0.15, 0.5, 1.0]
    )


# rrup at which pgv is close to target
def find_rrup(pgv_target):
    rrup = 50.0
    while True:
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        pgv = compute_gmm(faultprop, siteprop, GMM.Br_13, "PGV")[0]
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
def auto_time2(vm_corners, srf_corners, ds_multiplier):
    """Calculates the sim duration from the bounds of the vm and srf
    :param vm_corners: A list of (lon, lat) tuples
    :param srf_corners: A [4n*2] numpy array of (lon, lat) pairs where n is the number of planes in the srf
    :param ds_multiplier: An integer"""
    # S wave arrival time is determined by the distance from the srf centroid to the furthest corner
    s_wave_arrival = (
            max([
                geo.ll_dist(*corner, (srf_corners[:, 0].max + srf_corners[:, 0].min) / 2,
                            (srf_corners[:, 1].max + srf_corners[:, 1].min) / 2)
                for corner in vm_corners
            ]) / 3.2
    )
    # Rrup is determined by the largest vm corner to nearest srf corner distance
    siteprop.Rrup = max([
        min([geo.ll_dist(*corner, *s_corner) for s_corner in srf_corners])
        for corner in vm_corners
    ])
    # magnitude is in faultprop
    ds = compute_gmm(faultprop, siteprop, GMM.AS_16, "Ds595")[0]
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


def rrup2xylen(rrup, hh, points, rot=0, wd="."):
    """
    rrup: in km
    """
    lon_mid, lat_mid, dx_km, dy_km = gmt.region_fit_oblique(points, 90 - rot, wd=wd)
    # extend by wanted rrup
    min_x = geo.ll_shift(lat_mid, lon_mid, rrup + dx_km, 270 - rot)[::-1]
    max_x = geo.ll_shift(lat_mid, lon_mid, rrup + dx_km, 90 - rot)[::-1]
    min_y = geo.ll_shift(lat_mid, lon_mid, rrup + dy_km, 180 - rot)[::-1]
    max_y = geo.ll_shift(lat_mid, lon_mid, rrup + dy_km, 0 - rot)[::-1]
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
            }, "{}.yaml".format(vm_params)
        )


# get outer corners of a domain
def build_corners(origin, rot, xlen, ylen):
    # wanted xlen, ylen is at corners
    # approach answer at infinity / "I don't know the formula" algorithm

    # amount to shift from middle
    x_shift = xlen / 2.0
    y_shift = ylen / 2.0
    # x lengths will always be correct
    target_y = ylen

    # adjust middle ylen such that edge (corner to corner) ylen is correct
    while True:
        # upper and lower mid points
        t_u = geo.ll_shift(origin[1], origin[0], y_shift, rot)[::-1]
        t_l = geo.ll_shift(origin[1], origin[0], y_shift, rot + 180)[::-1]
        # bearings have changed at the top and bottom
        top_bearing = (geo.ll_bearing(t_u[0], t_u[1], origin[0], origin[1]) + 180) % 360
        bottom_bearing = geo.ll_bearing(t_l[0], t_l[1], origin[0], origin[1])
        # extend top and bottom middle line to make right edge
        c1 = geo.ll_shift(t_u[1], t_u[0], x_shift, top_bearing + 90)[::-1]
        c4 = geo.ll_shift(t_l[1], t_l[0], x_shift, bottom_bearing + 90)[::-1]
        # check right edge distance
        current_y = geo.ll_dist(c1[0], c1[1], c4[0], c4[1])
        if round(target_y, 10) != round(current_y, 10):
            y_shift += (target_y - current_y) / 2.0
        else:
            break
    # complete for left edge
    c2 = geo.ll_shift(t_u[1], t_u[0], x_shift, top_bearing - 90)[::-1]
    c3 = geo.ll_shift(t_l[1], t_l[0], x_shift, bottom_bearing - 90)[::-1]

    # at this point we have a perfect square (by corner distance)
    # c1 -> c4 == c2 -> c3 (right == left), c1 -> c2 == c3 -> c4 (top == bottom)
    return c1, c2, c3, c4


# maximum width along lines a and b
def reduce_domain(a0, a1, b0, b1, hh, space_srf, space_land, wd):
    # scan line paths follow lines a, b
    len_ab = geo.ll_dist(a0[0], a0[1], a1[0], a1[1])
    bearing_a = geo.ll_bearing(a0[0], a0[1], a1[0], a1[1])
    bearing_b = geo.ll_bearing(b0[0], b0[1], b1[0], b1[1])

    # determine rest of model parameters
    origin = geo.ll_mid(a0[0], a0[1], b1[0], b1[1])
    centre_top = geo.ll_mid(a1[0], a1[1], b1[0], b1[1])
    rot = geo.ll_bearing(origin[0], origin[1], centre_top[0], centre_top[1])
    xlen = geo.ll_dist(a0[0], a0[1], b0[0], b0[1])

    # number of scan lines accross domain (inclusive of edges)
    scanlines = 81
    # first scan, reduce east and west
    over_e = float("-inf") #None
    over_w = float("-inf") #None
    # store scan data for 2nd level processing
    scan_extremes = np.zeros((scanlines, 2, 2))
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint=True)):
        a = geo.ll_shift(a0[1], a0[0], x * len_ab, bearing_a)[::-1]
        b = geo.ll_shift(b0[1], b0[0], x * len_ab, bearing_b)[::-1]
        m = geo.ll_mid(a[0], a[1], b[0], b[1])
        # extend path a -> b, make sure we reach land for intersections
        bearing_p = geo.ll_bearing(m[0], m[1], a[0], a[1])
        ap = geo.ll_shift(m[1], m[0], 600, bearing_p)[::-1]
        bp = geo.ll_shift(m[1], m[0], 600, bearing_p + 180)[::-1]
        # mid to east edge, west edge
        m2ee = geo.ll_dist(m[0], m[1], a[0], a[1])
        m2we = geo.ll_dist(m[0], m[1], b[0], b[1])
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
        as_east = []
        as_west = []

        for c in range(len(isections)):
            if "%s/srf.path" % (wd) in icomps[c]:
                diff = space_srf
            elif NZ_LAND_OUTLINE in icomps[c]:
                diff = space_land
            # these bearings will be slightly wrong but don't need to be exact
            # TODO: adjust to follow same path as ap -> bp
            as_east.append(
                geo.ll_shift(isections[c][1], isections[c][0], diff, bearing_p)[::-1]
            )
            as_west.append(
                geo.ll_shift(isections[c][1], isections[c][0], diff, bearing_p + 180)[
                    ::-1
                ]
            )
        # find final extremes
        east = as_east[np.argmax(as_east, axis=0)[0]]
        west = as_west[np.argmin(as_west, axis=0)[0]]
        # for north/south (later on), use direct/un-shifted points
        scan_extremes[i] = (
            isections[np.argmax(isections, axis=0)[0]],
            isections[np.argmin(isections, axis=0)[0]],
        )

        # distance beyond east is negative
        m2e = geo.ll_dist(m[0], m[1], east[0], east[1])
        be = geo.ll_bearing(m[0], m[1], east[0], east[1])
        if abs(be - (bearing_p) % 360) > 90:
            m2e *= -1
        over_e = max(over_e, math.ceil((m2e - m2ee) / hh) * hh)
        # distance beyond west is negative
        m2w = geo.ll_dist(m[0], m[1], west[0], west[1])
        bw = geo.ll_bearing(m[0], m[1], west[0], west[1])
        if abs(bw - (bearing_p + 180) % 360) > 90:
            m2w *= -1
        over_w = max(over_w, math.ceil((m2w - m2we) / hh) * hh)

    # bring in sides
    xlen2 = xlen + over_e * (over_e < 0) + over_w * (over_w < 0)
    origin2 = geo.ll_shift(
        origin[1], origin[0], xlen2 / 2.0 - over_w * (over_w < 0) - xlen / 2.0, rot + 90
    )[::-1]
    a1, b1, b0, a0 = build_corners(origin2, rot, xlen2, len_ab)

    # corners may have moved
    bearing_a = geo.ll_bearing(a0[0], a0[1], a1[0], a1[1])
    bearing_b = geo.ll_bearing(b0[0], b0[1], b1[0], b1[1])
    # second scan, reduce north and south
    over_s = 0
    over_n = len_ab
    land_discovered = False
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint=True)):
        a = geo.ll_shift(a0[1], a0[0], x * len_ab, bearing_a)[::-1]
        b = geo.ll_shift(b0[1], b0[0], x * len_ab, bearing_b)[::-1]
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
        if not isect_inside and (not isect_east or not isect_west):
            # TODO: invert above if statement gracefully
            pass
        else:
            if not land_discovered:
                # shift bottom edge up
                over_s = math.floor((x * len_ab - 15) / hh) * hh
            land_discovered = True
            over_n = math.floor((len_ab - x * len_ab - 15) / hh) * hh

    # bring in top and bottom edge
    len_ab2 = len_ab - over_s * (over_s > 0) - over_n * (over_n > 0)
    origin2 = geo.ll_shift(
        origin2[1],
        origin2[0],
        len_ab2 / 2.0 + over_s * (over_s > 0) - len_ab / 2.0,
        rot,
    )[::-1]
    a1, b1, b0, a0 = build_corners(origin2, rot, xlen2, len_ab2)

    return a0, a1, b0, b1


def gen_vm(args, srf_meta, vm_params_dict, mag, ptemp):
    # store configs
    vm_dir = os.path.join(args.out_dir, srf_meta["name"])
    print(vm_dir)
    vm_params_dict["vm_dir"] = vm_dir
    nzvm_cfg = os.path.join(ptemp, "nzvm.cfg")
    vm_params_path = os.path.join(ptemp, "vm_params")
    # NZVM won't run if folder exists
    if os.path.exists(vm_dir):
        rmtree(vm_dir)
    save_vm_config(
        nzvm_cfg=nzvm_cfg,
        vm_params=vm_params_path,
        vm_dir=vm_dir,
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
        # save important files
        os.makedirs(vm_dir)
        move(nzvm_cfg, vm_dir)
        move("%s.yaml" % (vm_params_path), vm_dir)
        # generate a corners like NZVM would have
        with open("%s/VeloModCorners.txt" % (vm_params_dict["vm_dir"]), "wb") as c:
            c.write("> VM corners (python generated)\n".encode())
            c.write(vm_params_dict["path_mod"].encode())
        return

    # NZVM won't find resources if WD is not NZVM dir, stdout not MPROC friendly
    with open(os.path.join(ptemp, "NZVM.out"), "w") as logfile:
        nzvm_exe = Popen(
            [NZVM_BIN, nzvm_cfg], cwd=os.path.dirname(NZVM_BIN), stdout=logfile
        )
        nzvm_exe.communicate()

    # fix up directory contents
    move(os.path.join(vm_dir, "Velocity_Model", "rho3dfile.d"), vm_dir)
    move(os.path.join(vm_dir, "Velocity_Model", "vp3dfile.p"), vm_dir)
    move(os.path.join(vm_dir, "Velocity_Model", "vs3dfile.s"), vm_dir)
    move(os.path.join(vm_dir, "Log", "VeloModCorners.txt"), vm_dir)
    rmtree(os.path.join(vm_dir, "Velocity_Model"))
    rmtree(os.path.join(vm_dir, "Log"))
    move(nzvm_cfg, vm_dir)
    move("%s.yaml" % (vm_params_path), vm_dir)
    # create model_coords, model_bounds etc...
    gen_coords(vm_dir=vm_dir)
    # validate
    success, message = validate_vm(vm_dir)
    if success:
        sys.stderr.write("VM check OK: %s\n" % (vm_dir))
    else:
        sys.stderr.write("VM check BAD: %s\n" % (message))


def plot_vm(vm_params, srf_corners, mag, ptemp):
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
        p.points(
            "%s/VeloModCorners.txt" % (vm_params["vm_dir"]),
            fill="red",
            line=None,
            shape="c",
            size=0.05,
        )
    else:
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
    if vm_exists:
        p.png(dpi=200, clip=True, background="white", out_dir=vm_params["vm_dir"])
    else:
        p.png(
            dpi=200,
            clip=True,
            background="white",
            out_name=os.path.abspath(os.path.join(ptemp, os.pardir, vm_params["name"])),
        )


def create_vm(args, srf_meta):
    # temp directory for current process
    ptemp = mkdtemp(prefix="_tmp_%s_" % (srf_meta["name"]), dir=args.out_dir)

    # properties stored in classes (fault of external code)
    faultprop.Mw = srf_meta["mag"]
    faultprop.rake = srf_meta["rake"]
    faultprop.dip = srf_meta["dip"]
    # rrup to reach wanted PGV
    if args.pgv == -1.0:
        pgv = mag2pgv(faultprop.Mw)
    else:
        pgv = args.pgv
    rrup, pgv_actual = find_rrup(pgv)

    # original, unrotated vm
    bearing = 0
    origin, xlen0, ylen0 = rrup2xylen(
        rrup, args.hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=ptemp
    )
    o1, o2, o3, o4 = build_corners(origin, bearing, xlen0, ylen0)
    vm0_region = corners2region(o1, o2, o3, o4)
    plot_region = (
        vm0_region[0] - 1,
        vm0_region[1] + 1,
        vm0_region[2] - 1,
        vm0_region[3] + 1,
    )

    # proportion in ocean
    if args.no_optimise:
        land0 = 100
    else:
        land0 = vm_land(o1, o2, o3, o4, wd=ptemp)

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
        print("modifying %s" % (srf_meta["name"]))

        # rotation based on centre line bearing at this point
        l1 = centre_lon(vm0_region[2])
        l2 = centre_lon(vm0_region[3])
        mid = geo.ll_mid(l1, vm0_region[2], l2, vm0_region[3])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, vm0_region[3]))

        # wanted distance is at corners, not middle top to bottom
        # approach answer at infinity / "I don't know the formula" algorithm
        xlen1, ylen1 = rrup2xylen(
            rrup, args.hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=ptemp
        )[1:]
        try:
            c1, c2, c3, c4 = build_corners(origin, bearing, ylen1, xlen1)
        except ValueError:
            raise ValueError("Error for vm {}".format(srf_meta["name"]))

        # cut down ocean areas
        c4, c1, c3, c2 = reduce_domain(
            c4, c1, c3, c2, args.hh, args.space_srf, args.space_land, ptemp
        )
        origin = geo.ll_mid(c4[0], c4[1], c2[0], c2[1])
        ylen1 = math.ceil(geo.ll_dist(c4[0], c4[1], c1[0], c1[1]) / args.hh) * args.hh
        xlen1 = math.ceil(geo.ll_dist(c4[0], c4[1], c3[0], c3[1]) / args.hh) * args.hh

        # proportion in ocean
        land1 = vm_land(c1, c2, c3, c4, wd=ptemp)

        # adjust region to fit new corners
        plot_region = (
            min(c4[0], c1[0], c3[0], c2[0], plot_region[0]),
            max(c4[0], c1[0], c3[0], c2[0], plot_region[1]),
            min(c4[1], c1[1], c3[1], c2[1], plot_region[2]),
            max(c4[1], c1[1], c3[1], c2[1], plot_region[3]),
        )
    else:
        print("not modifying %s" % (srf_meta["name"]))
        # store original variables as final versions
        ylen1 = ylen0
        xlen1 = xlen0
        land1 = land0
        c1, c2, c3, c4 = o1, o2, o3, o4
    # not enough land in final domain
    if math.floor(land1) == 0:
        xlen1 = 0
        ylen1 = 0

    # zlen is independent from xlen and ylen
    zlen = round(auto_z(faultprop.Mw, srf_meta["dbottom"]) / args.hh) * args.hh
    # modified sim time
    vm_corners = (c1, c2, c3, c4)
    initial_time = auto_time2(vm_corners, np.concatenate(srf_meta["corners"], axis=0), 1.2)
    sim_time1 = (initial_time // args.dt) * args.dt

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
        gen_vm(args, srf_meta, vm_params, faultprop.Mw, ptemp)
    # plot results
    plot_vm(vm_params, srf_meta["corners"], faultprop.Mw, ptemp)

    # working dir cleanup, return info about VM
    rmtree(ptemp)
    return vm_params


def load_msgs(args):
    # returns list of appropriate srf metadata
    msgs = []
    faults = set()

    # add task for every info file
    info_files = glob(args.info_glob)
    for info in info_files:
        if not os.path.exists(info):
            print("SRF info file not found: %s" % (info))
            continue

        # name is unique and based on basename
        name = os.path.splitext(os.path.basename(info))[0].split("_")[0]

        if name in faults:
            continue
        faults.add(name)

        with h5open(info) as h:
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
                )
            )
    return msgs


def load_msgs_nhm(args):
    msgs = []

    with open(args.nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
    n_i = 0
    while n_i < len(nhm):
        n = nhm[n_i].strip()
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
            )
        )

    return msgs


def store_nhm_selection(out_dir, reports):
    selected = os.path.join(out_dir, "nhm_selection.txt")
    excluded = os.path.join(out_dir, "excluded.txt")
    with open(selected, "w") as sf:
        with open(excluded, "w") as ef:
            for r in reports:
                if r["xlen_mod"] == 0 or r["ylen_mod"] == 0 or r["zlen"] == 0:
                    ef.write("%s\n" % (r["name"]))
                else:
                    sf.write(
                        "%s %dr\n" % (r["name"], min(max(r["mag"] * 20 - 110, 10), 50))
                    )


def store_summary(table, info_store):
    # initialise table file
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


def load_args():
    parser = ArgumentParser()
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
        "-o", "--out-dir", help="directory to place outputs", default="autovm"
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
    arg("--novm", help="only generate parameters", action="store_true")
    arg("--vm-version", help="velocity model version to generate", default="1.65")
    arg(
        "--vm-topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )
    arg("--selection", help="also generate NHM selection file", action="store_true")
    arg("--no_optimise", help="Don't try and optimise the vm if it is off shore. Removes dependency on having GMT coastline data", action="store_true")
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)
    if not args.novm:
        if NZVM_BIN is None:
            sys.exit("""NZVM binary not in PATH
You can compile it from here: https://github.com/ucgmsim/Velocity-Model
Then add it to the PATH by either:
sudo cp /location/to/Velocity-Model/NZVM /usr/bin/
or:
export PATH=$PATH:/location/to/Velocity-Model""")

    return args


if __name__ == "__main__":
    args = load_args()

    # load wanted fault information
    if args.info_glob == "NHM":
        msg_list = load_msgs_nhm(args)
    else:
        msg_list = load_msgs(args)
    if len(msg_list) == 0:
        sys.exit("Found nothing to do.")

    # prepare to run
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    # proxy function for mapping
    def create_vm_star(args_meta):
        return create_vm(*args_meta)

    # distribute work
    if args.nproc > 1:
        p = Pool(processes=args.nproc)
        reports = p.map(create_vm_star, msg_list)
    else:
        # debug friendly alternative
        reports = [create_vm_star(msg) for msg in msg_list]

    # nhm selection formatted file and list of excluded VMs
    if args.selection:
        store_nhm_selection(args.out_dir, reports)
    # store summary
    store_summary(os.path.join(args.out_dir, "vminfo.csv"), reports)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(('maui', 'mahuika', 'wb', 'ni')): # Checks if is on the HPCF
        permission_cmd = ['chmod', 'g+rwXs', '-R', args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ['chgrp', constants.DEFAULT_ACCOUNT, '-R', args.out_dir]
        subprocess.call(group_cmd)

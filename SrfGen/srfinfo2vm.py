#!/usr/bin/env python
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example (set out_dir with -o or --out-dir):
mpirun -n 3 srfinfo2vm.py "Srf/*.info"

HELP:
./srfinfo2vm -h
"""
import argparse
from distutils.spawn import find_executable
from glob import glob
from typing import Dict

import math
import os
import platform
import subprocess
from shutil import rmtree, move
import sys
from tempfile import TemporaryDirectory, mk

from h5py import File as h5open
import numpy as np

# qcore library should already be in path
from qcore.simulation_structure import get_fault_from_realisation

from qcore import constants
from qcore import geo
from qcore import gmt
from qcore.validate_vm import validate_vm
from qcore.utils import dump_yaml

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

from SrfGen.createSRF import leonard, mom2mag, mag2mom
from SrfGen.gen_coords import gen_coords

script_dir = os.path.dirname(os.path.abspath(__file__))
NZ_CENTRE_LINE = os.path.join(script_dir, "NHM/res/centre.txt")
NZ_LAND_OUTLINE = os.path.join(script_dir, "NHM/res/rough_land.txt")
NZVM_BIN = find_executable("NZVM")


# default scaling relationship
def mag2pgv(mag):
    return np.interp(
        mag, [3.5, 4.1, 4.7, 5.2, 6.4, 7.5], [0.015, 0.0375, 0.075, 0.15, 0.5, 1.0]
    )


# rrup at which pgv is close to target
def find_rrup(pgv_target, faultprop, siteprop):
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
def auto_time2(xlen, ylen, ds_multiplier, faultprop, siteprop):
    rrup = max(xlen, ylen) / 2.0
    s_wave_arrival = rrup / 3.2
    siteprop.Rrup = rrup
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

def save_nzvm_config(
    nzvm_cfg=None,
    vm_dir=None,
    origin=(170, -40),
    rot=0,
    xlen=100,
    ylen=100,
    zmax=40,
    zmin=0,
    hh=0.4,
    min_vs=0.5,
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
    assert vm_dir is not None and not os.path.exists(vm_dir)
    with open(nzvm_cfg, "w") as vmd:
        vmd.write(
            "\n".join(
                [
                    "CALL_TYPE=GENERATE_VELOCITY_MOD",
                    "MODEL_VERSION={}".format(model_version),
                    "OUTPUT_DIR={}".format(vm_dir),
                    "ORIGIN_LAT={}".format(origin[1]),
                    "ORIGIN_LON={}".format(origin[0]),
                    "ORIGIN_ROT={}".format(rot),
                    "EXTENT_X={}".format(xlen),
                    "EXTENT_Y={}".format(ylen),
                    "EXTENT_ZMAX={}".format(zmax),
                    "EXTENT_ZMIN={}".format(zmin),
                    "EXTENT_Z_SPACING={}".format(hh),
                    "EXTENT_LATLON_SPACING={}".format(hh),
                    "MIN_VS={}".format(min_vs),
                    "TOPO_TYPE={}\n".format(topo_type),
                ]
            )
        )


def generate_vm_params(
    vm_gen_params,
    vm_dir,
    zmin=0,
    hh=0.4,
    min_vs=0.5,
    centroid_depth=7,
    code="rt",
    model_version="1.65",
    topo_type="BULLDOZED",
):
    model_lon, model_lat = vm_gen_params["origin"]
    xlen_mod = vm_gen_params["xlen_mod"]
    ylen_mod = vm_gen_params["ylen_mod"]
    zlen_mod = vm_gen_params["zlen_mod"]
    return {
        "mag": float(vm_gen_params["mag"]),
        "centroidDepth": float(centroid_depth),
        "MODEL_LAT": float(model_lat),
        "MODEL_LON": float(model_lon),
        "MODEL_ROT": float(vm_gen_params["bearing"]),
        "hh": hh,
        "min_vs": float(min_vs),
        "model_version": model_version,
        "topo_type": topo_type,
        "output_directory": os.path.basename(vm_dir),
        "code": code,
        "extent_x": float(xlen_mod),
        "extent_y": float(ylen_mod),
        "extent_zmax": float(zlen_mod),
        "extent_zmin": float(zmin),
        "sim_duration": float(vm_gen_params["sim_time_mod"]),
        "flo": min_vs / (5.0 * hh),
        "nx": int(round(float(xlen_mod) / hh)),
        "ny": int(round(float(ylen_mod) / hh)),
        "nz": int(round(float(zlen_mod - zmin) / hh)),
        "sufx": "_{}01-h{:.3f}".format(code, hh),
    }


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

    min_lon = min(a0[0], a1[0], b0[0], b1[0])
    max_lon = max(a0[0], a1[0], b0[0], b1[0])
    min_lat = min(a0[1], a1[1], b0[1], b1[1])
    max_lat = max(a0[1], a1[1], b0[1], b1[1])

    if max_lon - min_lon > 180 - 165 or max_lat -min_lat > (-33) - (-48):
        raise ValueError(
            "VM does not fit within NZ DEM bounds. {}, {}, {}, {}, rot: {}, origin {}".format(
                a0, a1, b0, b1, rot, origin
            )
        )

    diff_lat, diff_lon = 0, 0

    if min_lon < 165:
        diff_lon = 1.001*(165 - min_lon)
    elif max_lon > 180:
        diff_lon = 1.001*(180 - max_lon)

    if min_lat < -48:
        diff_lat = 1.001*(-48 - min_lat)
    elif max_lat > -33:
        diff_lat = 1.001*(-33 - max_lat)

    if diff_lat or diff_lon:
        a0 = (a0[0] + diff_lon, a0[1] + diff_lat)
        a1 = (a1[0] + diff_lon, a1[1] + diff_lat)
        b0 = (b0[0] + diff_lon, b0[1] + diff_lat)
        b1 = (b1[0] + diff_lon, b1[1] + diff_lat)
        
    return a0, a1, b0, b1


def gen_vm(min_vs, vm_topo, vm_version, novm, hh, out_dir, srf_meta, vm_gen_params, ptemp):
    # store configs
    vm_dir = os.path.join(out_dir, srf_meta["name"])
    nzvm_cfg = os.path.join(ptemp, "nzvm.cfg")
    vm_params_yaml_path = os.path.join(ptemp, "vm_params.yaml")
    # NZVM won't run if folder exists
    rmtree(vm_dir, ignore_errors=True)

    vm_params = generate_vm_params(
        vm_gen_params,
        vm_dir,
        hh=hh,
        min_vs=min_vs,
        centroid_depth=srf_meta["hdepth"],
        topo_type=vm_topo,
        model_version=vm_version,
    )
    dump_yaml(vm_params, vm_params_yaml_path)
    save_nzvm_config(
        nzvm_cfg=nzvm_cfg,
        vm_dir=vm_dir,
        origin=vm_gen_params["origin"],
        rot=vm_gen_params["bearing"],
        xlen=vm_gen_params["xlen_mod"],
        ylen=vm_gen_params["ylen_mod"],
        zmax=vm_gen_params["zlen_mod"],
        hh=hh,
        min_vs=min_vs,
        topo_type=vm_topo,
        model_version=vm_version,
    )

    if novm:
        # save important files
        os.makedirs(vm_dir)
        move(nzvm_cfg, vm_dir)
        move(vm_params_yaml_path, vm_dir)
        # generate a corners like NZVM would have
        with open("{}/VeloModCorners.txt".format(vm_dir), "wb") as c:
            c.write("> VM corners (python generated)\n".encode())
            c.write(vm_gen_params["path_mod"].encode())
    else:
        # NZVM won't find resources if WD is not NZVM dir, stdout not MPROC friendly
        with open(os.path.join(ptemp, "NZVM.out"), "w") as logfile:
            nzvm_exe = subprocess.Popen(
                [NZVM_BIN, nzvm_cfg], cwd=os.path.dirname(NZVM_BIN), stdout=logfile
            )
            nzvm_exe.communicate()

        # fix up directory contents
        files_to_move = [
            os.path.join(vm_dir, "Velocity_Model", "rho3dfile.d"),
            os.path.join(vm_dir, "Velocity_Model", "vp3dfile.p"),
            os.path.join(vm_dir, "Velocity_Model", "vs3dfile.s"),
            os.path.join(vm_dir, "Log", "VeloModCorners.txt"),
            nzvm_cfg,
            vm_params_yaml_path,
        ]
        for file_name in files_to_move:
            move(file_name, vm_dir)

        rmtree(os.path.join(vm_dir, "Velocity_Model"))
        rmtree(os.path.join(vm_dir, "Log"))

        # create model_coords, model_bounds etc...
        gen_coords(vm_dir=vm_dir)
        # validate
        success, message = validate_vm(vm_dir)
        if success:
            sys.stderr.write("VM check OK: {}\n".format(vm_dir))
        else:
            sys.stderr.write("VM check BAD: {}\n".format(message))


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


def create_vm_gen_params(temp_dir, srf_meta, space_land=5.0, space_srf=15.0, dt=0.005, hh=0.4, pgv=5.0):
    # temp directory for current process
    faultprop = Fault()
    siteprop = Site()

    # properties stored in classes (fault of external code)
    faultprop.Mw = srf_meta["mag"]
    faultprop.dip = srf_meta["dip"]
    faultprop.ztor = 0.0
    siteprop.vs30 = 500
    # rrup to reach wanted PGV
    if pgv == -1.0:
        pgv = mag2pgv(faultprop.Mw)
    rrup, pgv_actual = find_rrup(pgv, faultprop, siteprop)

    # original, unrotated vm
    bearing = 0
    origin, xlen0, ylen0 = rrup2xylen(
        rrup, hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=temp_dir
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
    land0 = vm_land(o1, o2, o3, o4, wd=temp_dir)

    # for plotting and calculating VM domain distance
    with open("{}/srf.path".format(temp_dir), "wb") as sp:
        for plane in srf_meta["corners"]:
            sp.write("> srf plane\n".encode())
            np.savetxt(sp, plane, fmt="%f")
            sp.write("%f %f\n".encode() % (tuple(plane[0])))

    # modify VM if necessary
    adjusted = False
    if faultprop.Mw >= 3.5 and land0 < 99:
        adjusted = True
        print("modifying {}".format(srf_meta["name"]))

        # rotation based on centre line bearing at this point
        l1 = centre_lon(vm0_region[2])
        l2 = centre_lon(vm0_region[3])
        mid = geo.ll_mid(l1, vm0_region[2], l2, vm0_region[3])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, vm0_region[3]))

        # wanted distance is at corners, not middle top to bottom
        # approach answer at infinity / "I don't know the formula" algorithm
        xlen1, ylen1 = rrup2xylen(
            rrup, hh, srf_meta["corners"].reshape((-1, 2)), rot=bearing, wd=temp_dir
        )[1:]
        try:
            c1, c2, c3, c4 = build_corners(origin, bearing, ylen1, xlen1)
        except ValueError:
            raise ValueError("Error for vm {}".format(srf_meta["name"]))

        # cut down ocean areas
        c4, c1, c3, c2 = reduce_domain(
            c4, c1, c3, c2, hh, space_srf, space_land, temp_dir
        )
        origin = geo.ll_mid(c4[0], c4[1], c2[0], c2[1])
        ylen1 = math.ceil(geo.ll_dist(c4[0], c4[1], c1[0], c1[1]) / hh) * hh
        xlen1 = math.ceil(geo.ll_dist(c4[0], c4[1], c3[0], c3[1]) / hh) * hh

        # proportion in ocean
        land1 = vm_land(c1, c2, c3, c4, wd=temp_dir)

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
    zlen = round(auto_z(faultprop.Mw, srf_meta["dbottom"]) / hh) * hh
    # modified sim time
    sim_time1 = (auto_time2(xlen1, ylen1, 1.2, faultprop, siteprop) // dt) * dt

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
        "path": "{} {}\n{} {}\n{} {}\n{} {}\n".format(
            o1[0], o1[1], o2[0], o2[1], o3[0], o3[1], o4[0], o4[1]
        ),
        "path_mod": "{} {}\n{} {}\n{} {}\n{} {}\n".format(
            c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]
        ),
    }


def store_nhm_selection(out_dir, vm_params):
    if min(vm_params["xlen_mod"], vm_params["ylen_mod"], vm_params["zlen_mod"]) > 0:
        selected = os.path.join(out_dir, "nhm_selection.txt")
        with open(selected, "a") as sf:
            sf.write(
                "{} {}r\n".format(vm_params["name"], min(max(vm_params["mag"] * 20 - 110, 10), 50))
            )
    else:
        excluded = os.path.join(out_dir, "excluded.txt")
        with open(excluded, "a") as ef:
            ef.write("{}\n".format(vm_params["name"]))


def store_summary(table, info_store):
    # initialise table file
    if not os.path.exists(table):
        with open(table, "w") as t:
            t.write(
                '"name","mw","plane depth (km, to bottom)",'
                '"vm depth (zlen, km)","xlen (km)",'
                '"ylen (km)","land (% cover)","adjusted vm depth (zlen, km)",'
                '"adjusted sim time (s)","adjusted xlen (km)",'
                '"adjusted ylen (km)","adjusted land (% cover)"\n'
            )
    with open(table, "a") as t:
        t.write(
            "{},{:.12g},{:g},{:g},{:g},{:g},{:.0f},{:g},{:g},{:g},{:g},{:.0f}\n".format(
                info_store["name"],
                info_store["mag"],
                info_store["dbottom"],
                info_store["zlen"],
                info_store["xlen"],
                info_store["ylen"],
                info_store["land"],
                info_store["zlen_mod"],
                info_store["sim_time_mod"],
                info_store["xlen_mod"],
                info_store["ylen_mod"],
                info_store["land_mod"],
            )
        )


def create_vm_from_srfinfo(out_dir, space_land, space_srf, dt, hh, pgv, min_vs, vm_topo, vm_version, novm, selection,
                           srf_meta):
    os.makedirs(out_dir, exist_ok=True)

    with TemporaryDirectory(prefix="_tmp_{}_".format(srf_meta["name"]), dir=out_dir) as temp_dir:
        vm_gen_params = create_vm_gen_params(
            temp_dir=temp_dir,
            srf_meta=srf_meta,
            space_land=space_land,
            space_srf=space_srf,
            dt=dt,
            hh=hh,
            pgv=pgv,
        )

        if min(vm_gen_params["xlen_mod"], vm_gen_params["ylen_mod"], vm_gen_params["zlen_mod"]) > 0:
            gen_vm(min_vs, vm_topo, vm_version, novm, hh, out_dir, srf_meta, vm_gen_params, temp_dir)

        plot_vm(vm_gen_params, srf_meta["corners"], vm_gen_params["mag"], temp_dir)

    # nhm selection formatted file and list of excluded VMs
    if selection:
        store_nhm_selection(out_dir, vm_gen_params)
    # store summary
    store_summary(os.path.join(out_dir, "vminfo.csv"), vm_gen_params)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(('maui', 'mahuika', 'wb', 'ni')):  # Checks if is on the HPCF
        permission_cmd = ['chmod', 'g+rwXs', '-R', out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ['chgrp', constants.DEFAULT_ACCOUNT, '-R', out_dir]
        subprocess.call(group_cmd)


def load_msg(info_file):
    # returns list of appropriate srf metadata

    # name is unique and based on basename
    name = os.path.splitext(os.path.basename(info_file))[0].split("_")[0]

    with h5open(info_file) as h:
        a = h.attrs
        try:
            mag = mom2mag(sum(map(mag2mom, a["mag"])))
        except TypeError:
            mag = a["mag"]
        return {
                    "name": name,
                    "dip": a["dip"][0],
                    "dbottom": a["dbottom"][0],
                    "corners": a["corners"],
                    "mag": mag,
                    "hdepth": a["hdepth"],
                    "rake": a["rake"],
                }


def load_msg_nhm(nhm_file: str, fault_name: str) -> Dict:

    with open(nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()

    n_i = 0
    while n_i < len(nhm) and nhm[n_i].strip() != fault_name:
        n_i += 13 + int(nhm[n_i + 11])
    if n_i >= len(nhm):
        return None
    name = nhm[n_i].strip()
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

    mag = leonard(rake, fwid * trace_length)
    return {
        "name": name,
        "dip": dip,
        "rake": rake,
        "dbottom": dbottom,
        "corners": corners,
        "mag": mag,
        "hdepth": dtop + 0.5 * (dbottom - dtop),
    }


def main():
    if NZVM_BIN is None:
        sys.exit("""NZVM binary not in PATH
            You can compile it from here: https://github.com/ucgmsim/Velocity-Model
            Then add it to the PATH by either:
            sudo cp /location/to/Velocity-Model/NZVM /usr/bin/
            or:
            export PATH=$PATH:/location/to/Velocity-Model""")

    parser = argparse.ArgumentParser()
    parser.add_argument("info_glob", help="info file selection expression. eg: Srf/*.info")
    parser.add_argument("out_dir", help="The VMs directory")
    parser.add_argument(
        "--nhm_file",
        help="path to NHM if using info_glob == 'NHM'",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "NHM", "NZ_FLTmodel_2010.txt"
        ),
    )
    parser.add_argument(
        "--pgv",
        help="max PGV at velocity model perimeter (estimated, cm/s)",
        type=float,
        default=5.0,
    )
    parser.add_argument("--hh", help="velocity model grid spacing (km)", type=float, default=0.4)
    parser.add_argument(
        "--dt",
        help="timestep to estimate simulation duration (s)",
        type=float,
        default=0.02,
    )
    parser.add_argument(
        "--space_land",
        help="min space between VM edge and land (km)",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--space_srf",
        help="min space between VM edge and SRF (km)",
        type=float,
        default=15.0,
    )
    parser.add_argument("--min_vs", help="for nzvm gen and flo (km/s)", type=float, default=0.5)
    parser.add_argument("--vm_version", help="velocity model version to generate", default="1.65")
    parser.add_argument(
        "--vm_topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )
    parser.add_argument("--selection", help="Append to the nhm selection file", action="store_true")
    parser.add_argument("--novm", help="only generate parameters", action="store_true")
    args = parser.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    # load wanted fault information

    info_files = glob(args.info_glob)

    use_nhm = args.info_glob == "NHM"

    faults = set()
    for info_file in info_files:

        name = get_fault_from_realisation(os.path.basename(info_file))

        if name in faults:
            continue
        faults.add(name)

        if use_nhm:
            srf_meta = load_msg_nhm(args.nhm_file, name)
            if srf_meta is None:
                raise ValueError("No NHM found for {}. Aborting velocity model generation".format(name))
        else:
            srf_meta = load_msg(info_file)

        create_vm_from_srfinfo(
            out_dir,
            args.space_land,
            args.space_srf,
            args.dt,
            args.hh,
            args.pgv,
            args.min_vs,
            args.vm_topo,
            args.vm_version,
            args.novm,
            args.selection,
            srf_meta,
        )


if __name__ == '__main__':
    main()

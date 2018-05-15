#!/usr/bin/env python2

import math
import os
from shutil import rmtree
from subprocess import Popen
import sys
from tempfile import mkdtemp
from time import time

import h5py as h5
from mpi4py import MPI
import numpy as np

# qcore library should already be in path
from qcore import geo
from qcore import gmt
# location containing Bradley_2010_Sa.py and AfshariStewart_2016_Ds.py
sys.path.append('/home/vap30/ucgmsim/post-processing/computations/')
from Bradley_2010_Sa import Bradley_2010_Sa
from AfshariStewart_2016_Ds import Afshari_Stewart_2016_Ds

NHM = 'NZ_FLTmodel_2010.txt'
# first entry index in NHM file
NHM_START = 15
PGV_TARGET = 5
NZ_CENTRE_LINE = 'res/centre.txt'
NZ_LAND_OUTLINE = 'res/rough_land.txt'
HH = 0.4
DT = 0.005
SPACE_LAND = 5
SPACE_SRF = 15
NZVM_BIN = '/home/vap30/ucgmsim/Velocity-Model/NZVM'
# call this as a script until the main function has a proper name
GEN_COORDS = '/nesi/projects/nesi00213/qcore/qcore/gen_coords.py'
# used to calculate flo (km / s)
MIN_VS = 0.5
# MPI - do not change
MASTER = 0

# Bradley_2010_Sa function requires passing parameters within classes
class siteprop:
    Rrup = None
    V30 = 500
    V30measured = 0
    Rx = Rrup
    Rjb = Rx
    Rtvz = 0
    # PGV
    period = -1
    Z1pt0 = math.exp(28.5 - (3.82 / 8.) * math.log(math.pow(V30, 8) \
            + math.pow(378.7, 8)))
    # for Ds
    defn = None
class faultprop:
    Mw = None
    Ztor = 0
    rake = None
    dip = None
    # for Ds
    rupture_type = None

# Leonard 2014 relationship between area and magnitude
def leonard(rake, A):
    # if dip slip else strike slip
    if round(rake % 360 / 90.) % 2:
        return 4.00 + math.log10(A)
    else:
        return 3.99 + math.log10(A)

# rrup at which pgv is close to target
def find_rrup():
    rrup = 50.0
    while True:
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        pgv = Bradley_2010_Sa(siteprop, faultprop)[0]
        # factor 0.02 is conservative step to avoid infinite looping
        if PGV_TARGET / pgv - 1 > 0.01:
            rrup -= rrup * 0.02
        elif PGV_TARGET / pgv - 1 < -0.01:
            rrup += rrup * 0.02
        else:
            break
    return rrup, pgv

# z extent based on magnitude, centroid depth
def auto_z(mag, depth):
    return round(10 + depth + (10 \
            * np.power((0.5 * np.power(10, (0.55 * mag - 1.2)) / depth), \
            (0.3))), 0)

# simulation time based on magnitude
def auto_time(mag):
    return np.power(10, max(1, 0.5 * mag - 1))

# simulation time based on area
def auto_time2(xlen, ylen, ds_multiplier):
    s_wave_arrival = max(xlen, ylen) / 3.2
    # alternative if rrup not available
    siteprop.Rrup = max(xlen / 2.0, ylen / 2.0)
    # magnitude is in faultprop
    ds = Afshari_Stewart_2016_Ds(siteprop, faultprop)[0]
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
        with h5.File(grd_file, 'r') as grd:
            return 100.0 * np.nansum(grd['z'][...]) / float(grd['z'].size)
    except IOError:
        print('INVALID GRD: %s' % (grd_file))
        # this should no longer happen when auto_dxy is used
        raise

# get longitude on NZ centre path at given latitude
def centre_lon(lat_target):
    # find closest points either side
    with open(NZ_CENTRE_LINE, 'r') as c:
        for line in c:
            ll = map(float, line.split())
            if ll[1] < lat_target:
                ll_bottom = ll
            else:
                ll_top = ll
                break
    try:
        # find rough proportion of distance
        lat_ratio = (lat_target - ll_bottom[1]) / (ll_top[1] - ll_bottom[1])
    except (UnboundLocalError, NameError):
        print('ERROR: VM edge out of lat bounds for centre line.')
        raise
    # consider same ratio along longitude difference to be correct
    return ll_bottom[0] + lat_ratio * (ll_top[0] - ll_bottom[0])

# return region that just fits 4 points
def corners2region(c1, c2, c3, c4):
    all_corners = np.array([c1, c2, c3, c4])
    x_min, y_min = np.min(all_corners, axis = 0)
    x_max, y_max = np.max(all_corners, axis = 0)
    return (x_min, x_max, y_min, y_max)

# create (closed) GMT path file from 4 points
def corners2path(c1, c2, c3, c4, output):
    with open(output, 'w') as o:
        o.write('>corners\n%s %s\n%s %s\n%s %s\n%s %s\n%s %s\n' \
                % (c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], \
                c4[0], c4[1], c1[0], c1[1]))

# get outer corners of a domain
def build_corners(origin, rot, xlen, ylen):
    # wanted xlen, ylen is at corners
    # approach answer at infinity / "I don't know the formula" algorithm

    # amount to shift from middle
    x_shift = xlen / 2.
    y_shift = ylen / 2.
    # x lengths will always be correct
    target_y = ylen

    # adjust middle ylen such that edge (corner to corner) ylen is correct
    while True:
        # upper and lower mid points
        t_u = geo.ll_shift(origin[1], origin[0], y_shift, rot)[::-1]
        t_l = geo.ll_shift(origin[1], origin[0], y_shift, rot + 180)[::-1]
        # bearings have changed at the top and bottom
        top_bearing = (geo.ll_bearing(t_u[0], t_u[1], origin[0], origin[1]) \
                + 180) % 360
        bottom_bearing = (geo.ll_bearing(t_l[0], t_l[1], origin[0], origin[1]))
        # extend top and bottom middle line to make right edge
        c1 = geo.ll_shift(t_u[1], t_u[0], x_shift, top_bearing + 90)[::-1]
        c4 = geo.ll_shift(t_l[1], t_l[0], x_shift, bottom_bearing + 90)[::-1]
        # check right edge distance
        current_y = geo.ll_dist(c1[0], c1[1], c4[0], c4[1])
        if round(target_y, 10) != round(current_y, 10):
            y_shift += (target_y - current_y) / 2.
        else:
            break
    # complete for left edge
    c2 = geo.ll_shift(t_u[1], t_u[0], x_shift, top_bearing - 90)[::-1]
    c3 = geo.ll_shift(t_l[1], t_l[0], x_shift, bottom_bearing - 90)[::-1]

    # at this point we have a perfect square (by corner distance)
    # c1 -> c4 == c2 -> c3 (right == left), c1 -> c2 == c3 -> c4 (top == bottom)
    return c1, c2, c3, c4

# maximum width along lines a and b
def reduce_domain(a0, a1, b0, b1, wd):
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
    over_e = None
    over_w = None
    # store scan data for 2nd level processing
    scan_extremes = np.zeros((scanlines, 2, 2))
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint = True)):
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
        geo.path_from_corners(corners = [ap, bp], \
                output = '%s/tempEXT.tmp' % (wd), \
                min_edge_points = 80, close = False)
        isections, icomps = gmt.intersections( \
                ['%s/srf.path' % (wd), NZ_LAND_OUTLINE, \
                '%s/tempEXT.tmp' % (wd)], items = True, \
                containing = '%s/tempEXT.tmp' % (wd))
        if len(isections) == 0:
            scan_extremes[i] = np.nan
            continue
        # shift different intersections by item-specific padding amount
        as_east = []
        as_west = []
        for c in xrange(len(isections)):
            if '%s/srf.path' % (wd) in icomps[c]:
                diff = SPACE_SRF
            elif NZ_LAND_OUTLINE in icomps[c]:
                diff = SPACE_LAND
            # these bearings will be slightly wrong but don't need to be exact
            # TODO: adjust to follow same path as ap -> bp
            as_east.append(geo.ll_shift(isections[c][1], isections[c][0], \
                    diff, bearing_p)[::-1])
            as_west.append(geo.ll_shift(isections[c][1], isections[c][0], \
                    diff, bearing_p + 180)[::-1])
        # find final extremes
        east = as_east[np.argmax(as_east, axis = 0)[0]]
        west = as_west[np.argmin(as_west, axis = 0)[0]]
        # for north/south (later on), use direct/un-shifted points
        scan_extremes[i] = isections[np.argmax(isections, axis = 0)[0]], \
                isections[np.argmin(isections, axis = 0)[0]]

        # distance beyond east is negative
        m2e = geo.ll_dist(m[0], m[1], east[0], east[1])
        be = geo.ll_bearing(m[0], m[1], east[0], east[1])
        if abs(be - (bearing_p) % 360) > 90:
            m2e *= -1
        over_e = max(over_e, math.ceil((m2e - m2ee) / HH) * HH)
        # distance beyond west is negative
        m2w = geo.ll_dist(m[0], m[1], west[0], west[1])
        bw = geo.ll_bearing(m[0], m[1], west[0], west[1])
        if abs(bw - (bearing_p + 180) % 360) > 90:
            m2w *= -1
        over_w = max(over_w, math.ceil((m2w - m2we) / HH) * HH)

    # bring in sides
    xlen2 = xlen + over_e * (over_e < 0) + over_w * (over_w < 0)
    origin2 = geo.ll_shift(origin[1], origin[0], \
            xlen2 / 2. - over_w * (over_w < 0) - xlen / 2., rot + 90)[::-1]
    a1, b1, b0, a0 = build_corners(origin2, rot, xlen2, len_ab)

    # corners may have moved
    bearing_a = geo.ll_bearing(a0[0], a0[1], a1[0], a1[1])
    bearing_b = geo.ll_bearing(b0[0], b0[1], b1[0], b1[1])
    # second scan, reduce north and south
    over_s = 0
    over_n = len_ab
    land_discovered = False
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint = True)):
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
                over_s = math.floor((x * len_ab - 15) / HH) * HH
            land_discovered = True
            over_n = math.floor((len_ab - x * len_ab - 15) / HH) * HH

    # bring in top and bottom edge
    len_ab2 = len_ab - over_s * (over_s > 0) - over_n * (over_n > 0)
    origin2 = geo.ll_shift(origin2[1], origin2[0], \
            len_ab2 / 2. + over_s * (over_s > 0) - len_ab / 2., rot)[::-1]
    a1, b1, b0, a0 = build_corners(origin2, rot, xlen2, len_ab2)
    return a0, a1, b0, b1

def create_vm(details):
    # temp directory for current process
    ptemp = mkdtemp(prefix = '_tmp_%s_' % (details['name']), \
            dir = details['out'])

    # properties stored in classes (fault of external code)
    faultprop.Mw = details['mag']
    faultprop.rake = details['rake']
    faultprop.dip = details['dip']
    # Karim: add 3km to depth if bottom >= 12km
    extend = details['dbottom'] >= 12
    # fault width (along dip)
    fwid = (details['dbottom'] - details['dtop'] + 3 * extend) \
        / math.sin(math.radians(faultprop.dip))
    # projected fault width (along dip direction)
    pwid = fwid * math.cos(math.radians(faultprop.dip))
    # Karim: fwid is different, update Mw
    if extend:
        trace_length = sum([geo.ll_dist(details['points'][plane][0], details['points'][plane][1], \
                details['points'][plane + 1][0], details['points'][plane + 1][1]) \
                for plane in xrange(details['n_pt'] - 1)])
        faultprop.Mw = leonard(faultprop.rake, fwid * trace_length)
    # rrup to reach wanted PGA
    rrup, pga_actual = find_rrup()

    # bottom edge of fault trace added after trace_length calculated
    details['points'].extend([geo.ll_shift(ll[1], ll[0], pwid, \
            details['dip_dir'])[::-1] for ll in details['points'][::-1]])
    # fault ll range (ll points at extremes)
    max_x, max_y = np.argmax(details['points'], axis = 0)
    min_x, min_y = np.argmin(details['points'], axis = 0)
    # extend by wanted rrup
    min_x = geo.ll_shift(details['points'][min_x][1], \
            details['points'][min_x][0], rrup, 270)[::-1]
    max_x = geo.ll_shift(details['points'][max_x][1], \
            details['points'][max_x][0], rrup, 90)[::-1]
    min_y = geo.ll_shift(details['points'][min_y][1], \
            details['points'][min_y][0], rrup, 180)[::-1]
    max_y = geo.ll_shift(details['points'][max_y][1], \
            details['points'][max_y][0], rrup, 0)[::-1]
    x_ext = math.ceil(geo.ll_dist(min_x[0], min_x[1], max_x[0], max_x[1]) \
            / HH) * HH
    y_ext = math.ceil(geo.ll_dist(min_y[0], min_y[1], max_y[0], max_y[1]) \
            / HH) * HH
    # vm region and plotting region
    ll_region0 = (min_x[0], max_x[0], min_y[1], max_y[1])
    ll_region = (min_x[0] - 1, max_x[0] + 1, min_y[1] - 1, max_y[1] + 1)
    # proper corners with equally sized edges
    bearing = 0
    mid0 = geo.ll_mid(min_x[0], min_y[1], max_x[0], max_y[1])
    o1, o2, o3, o4 = build_corners(mid0, bearing, x_ext, y_ext)

    # original domain has a ds_multiplier of 2 when calculating sim time
    sim_time = (auto_time2(x_ext, y_ext, 2) // DT) * DT

    # proportion in ocean for simple domain
    dxy = auto_dxy(corners2region(min_x, max_x, min_y, max_y))
    gmt.grd_mask('f', '%s/landmask.grd' % (ptemp), region = ll_region0, \
            dx = dxy, dy = dxy, wd = ptemp, outside = 0)
    land = grd_proportion('%s/landmask.grd' % (ptemp))

    # for plotting and calculating VM domain distance
    with open('%s/srf.path' % (ptemp), 'w') as sp:
        sp.write('\n'.join([' '.join(map(str, ll)) \
                for ll in details['points']]))
        sp.write('\n%s %s\n' \
                % (details['points'][0][0], details['points'][0][1]))

    # modify if necessary
    adjusted = False
    if faultprop.Mw >= 6.2 and land < 99:
        adjusted = True
        print('modifying %s' % (details['name']))

        # rotation based on centre line bearing at this point
        l1 = centre_lon(min_y[1])
        l2 = centre_lon(max_y[1])
        mid = geo.ll_mid(l1, min_y[1], l2, max_y[1])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, max_y[1]))

        # wanted distance is at corners, not middle top to bottom
        # approach answer at infinity / "I don't know the formula" algorithm
        m = geo.ll_mid(min_x[0], min_y[1], max_x[0], max_y[1])
        # TODO: recalculate xlen, ylen based on rrup and rotation
        c1, c2, c3, c4 = build_corners(m, bearing, y_ext, x_ext)

        # cut down ocean areas
        c4, c1, c3, c2 = reduce_domain(c4, c1, c3, c2, ptemp)

        # adjust region to fit new corners
        ll_region = (min(c4[0], c1[0], c3[0], c2[0], ll_region[0]), \
                max(c4[0], c1[0], c3[0], c2[0], ll_region[1]), \
                min(c4[1], c1[1], c3[1], c2[1], ll_region[2]), \
                max(c4[1], c1[1], c3[1], c2[1], ll_region[3]))

        # proportion in ocean
        ll_region1 = corners2region(c1, c2, c3, c4)
        dxy = auto_dxy(ll_region1)
        path_vm = '%s/landmask_vm.path' % (ptemp)
        grd_vm = '%s/landmask_vm.grd' % (ptemp)
        grd_land = '%s/landmask1.grd' % (ptemp)
        grd_and = '%s/landmask_AND.grd' % (ptemp)
        corners2path(c1, c2, c3, c4, path_vm)
        gmt.grd_mask('f', grd_land, \
                region = ll_region1, dx = dxy, dy = dxy, wd = ptemp)
        gmt.grd_mask(path_vm, grd_vm, \
                region = ll_region1, dx = dxy, dy = dxy, wd = ptemp)
        total_vm = grd_proportion(grd_vm)
        gmt.grdmath([grd_land, grd_vm, 'BITAND', '=', grd_and], wd = ptemp)
        try:
            land1 = grd_proportion(grd_and) / total_vm * 100
        except ZeroDivisionError:
            # this should not happen, bug in h5py or hdf5
            # auto_dxy should be a workaround
            land1 = 0

        mid0 = geo.ll_mid(c4[0], c4[1], c2[0], c2[1])
        ylen = math.ceil(geo.ll_dist(c4[0], c4[1], c1[0], c1[1]) / HH) * HH
        xlen = math.ceil(geo.ll_dist(c4[0], c4[1], c3[0], c3[1]) / HH) * HH
    else:
        print('not modifying %s' % (details['name']))
        ylen = y_ext
        xlen = x_ext
        land1 = land
    # not enough land in final domain
    if math.floor(land1) == 0:
        xlen = 0
        ylen = 0

    # modified sim time
    sim_time_mod = (auto_time2(xlen, ylen, 1) // DT) * DT
    # assuming hypocentre has been placed at 0.6 * depth
    centroid_depth = (details['dbottom'] - details['dtop']) * 0.6 \
            + details['dtop']
    zlen = round(auto_z(faultprop.Mw, centroid_depth) / HH) * HH

    if xlen == 0 or ylen == 0 or zlen == 0:
        # clean and only return metadata
        rmtree(ptemp)
        return {'name':details['name'], 'mag':faultprop.Mw, \
                'dbottom':details['dbottom'], \
                'zlen':zlen, 'sim_time':sim_time, \
                'xlen':x_ext, 'ylen':y_ext, 'land':land, \
                'zlen_mod':zlen, 'sim_time_mod':sim_time_mod, \
                'xlen_mod':xlen, 'ylen_mod':ylen, 'land_mod':land1}

    # store info
    ###
    ### XXX: params_vel depends on gen_coords.py
    ###
    vm_out = os.path.join(details['out'], details['name'])
    nzvm_cfg = os.path.join(ptemp, 'nzvm.cfg')
    params_vel = os.path.join(ptemp, 'params_vel.py')
    with open(nzvm_cfg, 'w') as vmd:
        vmd.write('\n'.join(['CALL_TYPE=GENERATE_VELOCITY_MOD', \
                'MODEL_VERSION=1.65', \
                'OUTPUT_DIR=%s' % (vm_out), \
                'ORIGIN_LAT=%s' % (mid0[1]), \
                'ORIGIN_LON=%s' % (mid0[0]), \
                'ORIGIN_ROT=%s' % (bearing), \
                'EXTENT_X=%s' % (xlen), \
                'EXTENT_Y=%s' % (ylen), \
                'EXTENT_ZMAX=%s' % (zlen), \
                'EXTENT_ZMIN=0.0', \
                'EXTENT_Z_SPACING=%s' % (HH), \
                'EXTENT_LATLON_SPACING=%s' % (HH), \
                'MIN_VS=0.5', \
                'TOPO_TYPE=BULLDOZED\n']))
    with open(params_vel, 'w') as pv:
        pv.write('\n'.join(['mag = "%s"' % (faultprop.Mw), \
                'centroidDepth = "%s"' % (centroid_depth), \
                'MODEL_LAT = "%s"' % (mid0[1]), \
                'MODEL_LON = "%s"' % (mid0[0]), \
                'MODEL_ROT = "%s"' % (bearing), \
                'hh = "%s"' % (HH), \
                'min_vs = "0.5"', \
                'model_version = "1.65"', \
                'topo_type = "BULLDOZED"', \
                'output_directory = "%s"' % (details['name']), \
                'extracted_slice_parameters_directory'
                ' = "SliceParametersNZ/SliceParametersExtracted.txt"', \
                'code = "rt"', \
                'extent_x = "%s"' % (xlen), \
                'extent_y = "%s"' % (ylen), \
                'extent_zmax = "%s"' % (zlen), \
                'extent_zmin = "0.0"', \
                'sim_duration = "%s"' % (sim_time * (not adjusted) \
                        + sim_time_mod * adjusted), \
                'flo = "%s"' % (MIN_VS / (5.0 * HH)), \
                'nx = "%s"' % (int(round(xlen / HH))), \
                'ny = "%s"' % (int(round(ylen / HH))), \
                'nz = "%s"' % (int(round(zlen / HH))), \
                'sufx = "_rt01-h%.3f"' % (HH)]))

    # create VM
    # limitation of NZVM code - won't run if folder exists
    if os.path.exists(vm_out):
        rmtree(vm_out)
    # limitation of NZVM code - won't find resources if wd is not NZVM dir
    # progress report is not multiprocess friendly (in stdout)
    with open('/dev/null', 'a') as sink:
        Popen([NZVM_BIN, nzvm_cfg], cwd = os.path.dirname(NZVM_BIN), \
                stdout = sink).wait()
    # fix up directory contents
    os.rename('%s/Velocity_Model/rho3dfile.d' % (vm_out), \
            '%s/rho3dfile.d' % (vm_out))
    os.rename('%s/Velocity_Model/vp3dfile.p' % (vm_out), \
            '%s/vp3dfile.p' % (vm_out))
    os.rename('%s/Velocity_Model/vs3dfile.s' % (vm_out), \
            '%s/vs3dfile.s' % (vm_out))
    os.rename('%s/Log/VeloModCorners.txt' % (vm_out), \
            '%s/VeloModCorners.txt' % (vm_out))
    rmtree(os.path.join(vm_out, 'Velocity_Model'))
    rmtree(os.path.join(vm_out, 'Log'))
    os.rename(nzvm_cfg, os.path.join(vm_out, os.path.basename(nzvm_cfg)))
    os.rename(params_vel, os.path.join(vm_out, os.path.basename(params_vel)))
    # create model_coords, model_bounds etc...
    # debug info unwanted (in stdout instead of stderr)
    with open('/dev/null', 'a') as sink:
        Popen([GEN_COORDS, vm_out], stdout = sink).wait()

    ###
    ### PLOT
    ###
    p = gmt.GMTPlot(os.path.join(vm_out, 'optimisation.ps'))
    p.spacial('M', ll_region, sizing = 7)
    p.coastlines()
    # filled slip area
    p.path('\n'.join([' '.join(map(str, ll)) for ll in details['points']]), \
            is_file = False, close = True, fill = 'yellow', split = '-')
    # top edge
    p.path('\n'.join([' '.join(map(str, ll)) \
            for ll in details['points'][:len(details['points']) / 2]]), is_file = False)
    # vm domain (simple and adjusted)
    p.path('%s %s\n%s %s\n%s %s\n%s %s' % (o1[0], o1[1], o2[0], o2[1], \
            o3[0], o3[1], o4[0], o4[1]), is_file = False, close = True, \
            fill = 'black@95')
    if adjusted:
        p.path('%s %s\n%s %s\n%s %s\n%s %s' % (c1[0], c1[1], c2[0], c2[1], \
                c3[0], c3[1], c4[0], c4[1]), is_file = False, close = True, \
                fill = 'black@95', split = '-', width = '1.0p')
    # info text for simple and adjusted domains
    p.text((ll_region[0] + ll_region[1]) / 2., ll_region[3], \
            'Mw: %.2f X: %.0fkm, Y: %.0fkm, land: %.0f%%' \
            % (faultprop.Mw, x_ext, y_ext, land), align = 'CT', dy = -0.1, \
            box_fill = 'white@50')
    if adjusted:
        p.text((ll_region[0] + ll_region[1]) / 2., ll_region[3], \
                'MODIFIED land: %.0f%%' \
                % (land1), align = 'CT', dy = -0.25, \
                box_fill = 'white@50')
    # for testing, show land mask cover
    # land outlines blue, nz centre line (for bearing calculation) red
    p.path(NZ_LAND_OUTLINE, is_file = True, close = False, \
            colour = 'blue', width = '0.2p')
    p.path(NZ_CENTRE_LINE, is_file = True, close = False, \
            colour = 'red', width = '0.2p')
    # actual corners retrieved from NZVM output
    p.points('%s/VeloModCorners.txt' % (vm_out), \
            fill = 'red', line = None, shape = 'c', size = 0.05)
    # store PNG, remove PS
    p.finalise()
    p.png(dpi = 200, clip = True, background = 'white')
    os.remove(p.pspath)

    # internal clean
    for tmp in ['gmt.conf', 'gmt.history']:
        if os.path.exists(os.path.join(vm_out, tmp)):
            os.remove(os.path.join(vm_out, tmp))
    # working dir cleanup
    rmtree(ptemp)

    # return vm info
    return {'name':details['name'], 'mag':faultprop.Mw, \
            'dbottom':details['dbottom'], \
            'zlen':zlen, 'sim_time':sim_time, \
            'xlen':x_ext, 'ylen':y_ext, 'land':land, \
            'zlen_mod':zlen, 'sim_time_mod':sim_time_mod, \
            'xlen_mod':xlen, 'ylen_mod':ylen, 'land_mod':land1}

def load_msgs(targets):
    msgs = []
    # PARAMETER 2: output folder (optional)
    if len(sys.argv) > 2:
        out = os.path.abspath(sys.argv[2])
    else:
        out = os.path.abspath('autovm')
    # make sure output directory exists
    if not os.path.exists(out):
        os.makedirs(out)

    # loop through faults
    with open(NHM, 'r') as nf:
        db = nf.readlines()
    dbi = NHM_START
    dbl = len(db)
    while dbi < dbl:
        name = db[dbi].strip()
        n_pt = int(db[dbi + 11])
        if targets != None and name not in targets:
            dbi += 13 + n_pt
            continue
        # retrieve required properties
        dip = float(db[dbi + 3].split()[0])
        dip_dir = float(db[dbi + 4])
        rake = float(db[dbi + 5])
        dbottom = float(db[dbi + 6].split()[0])
        dtop = float(db[dbi + 7].split()[0])
        mag = float(db[dbi + 10].split()[0])
        pts = [map(float, ll.split()) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
        # dictionary easy to retrieve, even with future changes
        msgs.append({'name':name, 'dip':dip, 'dip_dir':dip_dir, 'rake':rake, \
                'dbottom':dbottom, 'dtop':dtop, 'mag':mag, \
                'n_pt':n_pt, 'points':pts, 'out':out})
        # next fault
        dbi += 13 + n_pt
    return msgs

def store_summary(table, info_store):
    # initialise table file
    with open(table, 'w') as t:
        t.write('"name","mw","plane depth (km, to bottom)",'
                '"vm depth (zlen, km)","sim time (s)","xlen (km)",'
                '"ylen (km)","land (% cover)","adjusted vm depth (zlen, km)",'
                '"adjusted sim time (s)","adjusted xlen (km)",'
                '"adjusted ylen (km)","adjusted land (% cover)"\n')

        for p in info_store:
            for i in p:
                t.write('%s,%s,%s,%s,%s,%s,%s,%.0f,%s,%s,%s,%s,%.0f\n' \
                        % (i['name'], i['mag'], i['dbottom'], \
                        i['zlen'], i['sim_time'], \
                        i['xlen'], i['ylen'], i['land'], \
                        i['zlen_mod'], i['sim_time_mod'], \
                        i['xlen_mod'], i['ylen_mod'], i['land_mod']))
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == MASTER:
    print("master", rank)
###
### MASTER
###
    if len(sys.argv) > 1:
        # clock start
        t_start = MPI.Wtime()

        # PARAMETER 1: name of nhm selection file
        if sys.argv[1] == 'ALL':
            targets = None
        else:
            with open(sys.argv[1], 'r') as select:
                targets = [s[0] for s in map(str.split, select.readlines())]

        # load wanted fault information
        msg_list = load_msgs(targets)
        if len(msg_list) == 0:
            print('No matches found.')
            sys.exit(1)
        out = msg_list[0]['out']
        table = '%s/vminfo.csv' % (out)
        nproc = size - 1
        # # do not fully load a machine with many cores
        # nproc_max = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        # if nproc_max > 8:
        #     nproc_max = int(round(0.8 * nproc_max))
        # # not more processes than jobs
        # nproc = min(len(msg_list), nproc_max)
        # # spawn slaves
        # comm = MPI.COMM_WORLD.Spawn(
        #     sys.executable, args = [sys.argv[0]], maxprocs = nproc)
        status = MPI.Status()
        # distribute work to slaves who ask
        while nproc:
            # slave asking for work
            comm.recv(source = MPI.ANY_SOURCE, status = status)
            slave_id = status.Get_source()

            # no more jobs? kill signal
            if len(msg_list) == 0:
                msg_list.append(StopIteration)
                nproc -= 1

            # next job
            msg = msg_list[0]
            del(msg_list[0])
            comm.send(obj = msg, dest = slave_id)

        # gather, process reports from slaves
        reports = comm.gather([],root = MASTER)
        # print("reports",reports)
        store_summary(table, reports)
        # stop mpi
        # comm.Disconnect()

###
### SLAVE
###
else:
    print("slave", rank)
    # connect to parent
    # try:
    #     comm = MPI.Comm.Get_parent()
    #     rank = comm.Get_rank()
    # except MPI.Exception:
    #     print('First parameter is fault selection file.')
    #     print('Use "ALL" to loop through all faults.')
    #     print('First parameter not found.')
    #     print('Alternatively MPI cannot connect to parent.')
    #     sys.exit(1)

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):
        # task timer
        t0 = MPI.Wtime()
        # run
        # try:
        details = create_vm(task)
        # except:
        #     print('FAILED TASK: %s' % (task['name']))
        #     continue
        # store
        details['time'] = MPI.Wtime() - t0
        logbook.append(details)

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    # comm.Disconnect()

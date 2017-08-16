#!/usr/bin/env python2

import math
import os
import sys
from time import time

import h5py as h5
import numpy as np

import geo
import gmt
sys.path.append('/home/vap30/ucgmsim/post-processing/computations/')
from Bradley_2010_Sa import Bradley_2010_Sa
from AfshariStewart_2016_Ds import Afshari_Stewart_2016_Ds

nhm = 'NZ_FLTmodel_2010.txt'
pgv_target = 5
centre = 'res/centre.txt'
hh = 0.4
SPACE_LAND = 5
SPACE_SRF = 15

if len(sys.argv) > 1:
    target = sys.argv[1]
else:
    print('First parameter is name of fault in database.')
    print('Use "ALL" to loop through all faults.')
    sys.exit(1)
if target == 'ALL':
    target = None
if len(sys.argv) > 2:
    out = os.path.abspath(sys.argv[2])
else:
    out = os.path.abspath('autovm')
table = '%s/vminfo.csv' % (out)

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
    Z1pt0 = math.exp(28.5 - (3.82 / 8.) * math.log(math.pow(V30, 8) + math.pow(378.7, 8)))
    # for Ds
    defn = None
class faultprop:
    Mw = None
    Ztor = 0
    rake = None
    dip = None
    # for Ds
    rupture_type = None

# Leonard 2014 Relations
def leonard(rake, A):
    # if dip slip else strike slip
    if round(rake % 360 / 90.) % 2:
        return 4.19 + math.log10(A)
    else:
        return 4.18 + math.log10(A)

# rrup at which pgv is close to target
def find_rrup():
    rrup = 50.0
    while True:
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        pgv = Bradley_2010_Sa(siteprop, faultprop)[0]
        if pgv_target/pgv - 1 > 0.01:
            rrup -= rrup * 0.02
        elif pgv_target/pgv - 1 < -0.01:
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
# other variations
def auto_time2(xlen, ylen, ds_multiplier):
    s_wave_arrival = max(xlen, ylen) / 3.2
    # alternative if rrup not available
    siteprop.Rrup = max(xlen / 2.0, ylen / 2.0)
    # magnitude is in faultprop
    ds = Afshari_Stewart_2016_Ds(siteprop, faultprop)[0]
    return s_wave_arrival + ds_multiplier * ds

# keep dx and dy small enough relative to domain
def auto_dxy(region):
    x_diff = region[1] - region[0]
    y_diff = region[3] - region[2]

    # start off at roughly 1km
    dxy = 0.01
    # this formula may not be exact, should actually be more if anything (safe)
    while math.ceil(x_diff / dxy) * math.ceil(y_diff / dxy) < 16384:
        dxy *= 0.5
    return dxy

# proportion of gmt mask file filled (0 -> 1)
def grd_proportion(grd_file):
    try:
        with h5.File(grd_file, 'r') as grd:
            return 100.0 * np.nansum(grd['z'][...]) / float(grd['z'].size)
    except IOError:
        print('Warning, Invalid GRD')
        return 0.0

# get longitude on NZ centre path at latitude
def centre_lon(lat_target):
    # find closest points either side
    with open(centre, 'r') as c:
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
        print('WARNING: VM edge out of lat bounds for centre line.')
        raise
    # consider same ratio along longitude differenc to be correct
    return ll_bottom[0] + lat_ratio * (ll_top[0] - ll_bottom[0])

# 
def corners2region(c1, c2, c3, c4):
    all_corners = np.array([c1, c2, c3, c4])
    x_min, y_min = np.min(all_corners, axis = 0)
    x_max, y_max = np.max(all_corners, axis = 0)
    return (x_min, x_max, y_min, y_max)

#
def corners2path(c1, c2, c3, c4, output):
    with open(output, 'w') as o:
        o.write('>corners\n%s %s\n%s %s\n%s %s\n%s %s\n%s %s\n' \
                % (c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], \
                c4[0], c4[1], c1[0], c1[1]))

# maximum width along lines a and b
def max_width(a0, a1, b0, b1, m0, m1):
    len_a = geo.ll_dist(a0[0], a0[1], a1[0], a1[1])
    len_b = geo.ll_dist(b0[0], b0[1], b1[0], b1[1])
    len_m = geo.ll_dist(m0[0], m0[1], m1[0], m1[1])
    bearing_a = geo.ll_bearing(a0[0], a0[1], a1[0], a1[1])
    bearing_b = geo.ll_bearing(b0[0], b0[1], b1[0], b1[1])
    bearing_m = geo.ll_bearing(m0[0], m0[1], m1[0], m1[1])

    # number of scan lines accross domain (inclusive of edges)
    scanlines = 81
    # first scan, reduce east and west
    over_e = None
    over_w = None
    # store scan data for 2nd level processing
    scan_extremes = np.zeros((scanlines, 2, 2))
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint = True)):
        m = geo.ll_shift(m0[1], m0[0], x * len_m, bearing_m)[::-1]
        a = geo.ll_shift(a0[1], a0[0], x * len_a, bearing_a)[::-1]
        b = geo.ll_shift(b0[1], b0[0], x * len_b, bearing_b)[::-1]
        bearing_p = geo.ll_bearing(m[0], m[1], a[0], a[1])
        ap = geo.ll_shift(m[1], m[0], 600, bearing_p)[::-1]
        bp = geo.ll_shift(m[1], m[0], 600, bearing_p + 180)[::-1]
        # mid to east edge, west edge
        m2ee = geo.ll_dist(m[0], m[1], a[0], a[1])
        m2we = geo.ll_dist(m[0], m[1], b[0], b[1])
        # gmt spacial is not great-circle path connective
        geo.path_from_corners(corners = [ap, bp], \
                output = '%s/tempEXT.tmp' % (out), \
                min_edge_points = 80, close = False)
        isections, icomps = gmt.intersections( \
                ['%s/srf.path' % (out), 'res/rough_land.txt', \
                '%s/tempEXT.tmp' % (out)], items = True, \
                containing = '%s/tempEXT.tmp' % (out))
        if len(isections) == 0:
            scan_extremes[i] = np.nan
            continue
        # shift different intersections by item-specific padding amount
        as_east = []
        as_west = []
        for c in xrange(len(isections)):
            if '%s/srf.path' % (out) in icomps[c]:
                diff = SPACE_SRF
            elif 'res/rough_land.txt' in icomps[c]:
                diff = SPACE_LAND
            # these bearings will be slightly wrong but don't need to be exact
            as_east.append(geo.ll_shift(isections[c][1], isections[c][0], \
                    diff, bearing + 90)[::-1])
            as_west.append(geo.ll_shift(isections[c][1], isections[c][0], \
                    diff, bearing - 90)[::-1])
        # find final extremes
        east = as_east[np.argmax(as_east, axis = 0)[0]]
        west = as_west[np.argmin(as_west, axis = 0)[0]]
        # for north/south, use direct/un-shifted points
        scan_extremes[i] = isections[np.argmax(isections, axis = 0)[0]], \
                isections[np.argmin(isections, axis = 0)[0]]

        #
        m2e = geo.ll_dist(m[0], m[1], east[0], east[1])
        be = geo.ll_bearing(m[0], m[1], east[0], east[1])
        if abs(be - (bearing + 90) % 360) > 90:
            m2e *= -1
        over_e = max(over_e, m2e - m2ee)
        #
        m2w = geo.ll_dist(m[0], m[1], west[0], west[1])
        bw = geo.ll_bearing(m[0], m[1], west[0], west[1])
        if abs(bw - (bearing - 90) % 360) > 90:
            m2w *= -1
        over_w = max(over_w, m2w - m2we)

    # shift sides
    if over_e != None and over_e < 0:
        a0 = geo.ll_shift(a0[1], a0[0], over_e, bearing + 90)[::-1]
        a1 = geo.ll_shift(a1[1], a1[0], over_e, bearing + 90)[::-1]
    if over_e != None and over_w < 0:
        b0 = geo.ll_shift(b0[1], b0[0], over_w, bearing - 90)[::-1]
        b1 = geo.ll_shift(b1[1], b1[0], over_w, bearing - 90)[::-1]

    # corners have moved
    len_a = geo.ll_dist(a0[0], a0[1], a1[0], a1[1])
    len_b = geo.ll_dist(b0[0], b0[1], b1[0], b1[1])
    bearing_a = geo.ll_bearing(a0[0], a0[1], a1[0], a1[1])
    bearing_b = geo.ll_bearing(b0[0], b0[1], b1[0], b1[1])
    # second scan, reduce north and south
    over_s = None
    over_n = None
    land_discovered = False
    for i, x in enumerate(np.linspace(0, 1, scanlines, endpoint = True)):
        a = geo.ll_shift(a0[1], a0[0], x * len_a, bearing_a)[::-1]
        b = geo.ll_shift(b0[1], b0[0], x * len_b, bearing_b)[::-1]
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
            if not land_discovered:
                # shift bottom edge up
                over_s = x * len_m - 15
        else:
            land_discovered = True
            over_n = x * len_m + 15

    # shift top and bottom edge
    if over_n != None and over_n < len_a:
        a1 = geo.ll_shift(a0[1], a0[0], over_n, bearing_a)[::-1]
        b1 = geo.ll_shift(b0[1], b0[0], over_n, bearing_b)[::-1]
    if over_s != None and over_s > 0:
        a0 = geo.ll_shift(a0[1], a0[0], over_s, bearing_a)[::-1]
        b0 = geo.ll_shift(b0[1], b0[0], over_s, bearing_b)[::-1]
    return a0, a1, b0, b1

# make sure output directory exists for images
if not os.path.exists(out):
    os.makedirs(out)

# clear table file
with open(table, 'w') as t:
    t.write('"name","mw","plane depth (km, to bottom)","vm depth (zlen, km)","sim time (s)",'
            '"xlen (km)","ylen (km)","land (% cover)","adjusted sim time (s)",'
            '"adjusted xlen (km)","adjusted ylen (km)","adjusted land (% cover)"\n')

# loop through faults
with open(nhm, 'r') as nf:
    db = nf.readlines()
    dbi = 15
    dbl = len(db)
while dbi < dbl:
    name = db[dbi].strip()
    n_pt = int(db[dbi + 11])
    if target != None and name != target:
        dbi += 13 + n_pt
        continue
    faultprop.Mw = float(db[dbi + 10].split()[0])
    faultprop.rake = float(db[dbi + 5])
    faultprop.dip = float(db[dbi + 3].split()[0])
    rrup, pga_actual = find_rrup()

    # fault width (along dip), projected width (along dip direction)
    fwid = (float(db[dbi + 6].split()[0]) - float(db[dbi + 7].split()[0])) \
            / math.sin(math.radians(faultprop.dip))
    pwid = fwid * math.cos(math.radians(faultprop.dip))
    dip_dir = float(db[dbi + 4])

    # top of fault trace
    pts = [map(float, ll.split()) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
    # bottom of fault trace
    pts.extend([geo.ll_shift(ll[1], ll[0], pwid, dip_dir)[::-1] for ll in pts[::-1]])
    # fault ll range (ll points at extremes)
    max_x, max_y = np.argmax(pts, axis = 0)
    min_x, min_y = np.argmin(pts, axis = 0)
    # extend by wanted rrup
    min_x = geo.ll_shift(pts[min_x][1], pts[min_x][0], rrup, 270)[::-1]
    max_x = geo.ll_shift(pts[max_x][1], pts[max_x][0], rrup, 90)[::-1]
    min_y = geo.ll_shift(pts[min_y][1], pts[min_y][0], rrup, 180)[::-1]
    max_y = geo.ll_shift(pts[max_y][1], pts[max_y][0], rrup, 0)[::-1]
    x_ext = geo.ll_dist(min_x[0], min_x[1], max_x[0], max_x[1])
    y_ext = geo.ll_dist(min_y[0], min_y[1], max_y[0], max_y[1])
    # vm region and plotting region
    ll_region0 = (min_x[0], max_x[0], min_y[1], max_y[1])
    ll_region = (min_x[0] - 1, max_x[0] + 1, min_y[1] - 1, max_y[1] + 1)

    # other
    sim_time = auto_time(faultprop.Mw)
    sim_time = auto_time2(x_ext, y_ext, 2)

    # proportion in ocean
    dxy = auto_dxy(corners2region(min_x, max_x, min_y, max_y))
    gmt.grd_mask('f', '%s/landmask.grd' % (out), region = ll_region0, dx = dxy, dy = dxy, wd = out, outside = 0)
    land = grd_proportion('%s/landmask.grd' % (out))

    with open('%s/srf.path' % (out), 'w') as sp:
        sp.write('\n'.join([' '.join(map(str, ll)) for ll in pts]))
        sp.write('\n%s %s\n' % (pts[0][0], pts[0][1]))

    # modify if necessary
    adjusted = False
    if faultprop.Mw >= 6.2 and land < 90:
        adjusted = True
        print 'modifying', name
        # rotate based on centre line bearing
        l1 = centre_lon(min_y[1])
        l2 = centre_lon(max_y[1])
        mid = geo.ll_mid(l1, min_y[1], l2, max_y[1])
        mid0 = geo.ll_mid(min_x[0], min_y[1], max_x[0], max_y[1])
        bearing = round(geo.ll_bearing(mid[0], mid[1], l2, max_y[1]))
        t_u = geo.ll_shift(mid0[1], mid0[0], x_ext / 2., bearing)[::-1]
        t_l = geo.ll_shift(mid0[1], mid0[0], x_ext / 2., bearing + 180)[::-1]
        c1 = geo.ll_shift(t_u[1], t_u[0], y_ext / 2., bearing + 90)[::-1]
        c2 = geo.ll_shift(t_u[1], t_u[0], y_ext / 2., bearing - 90)[::-1]
        c3 = geo.ll_shift(t_l[1], t_l[0], y_ext / 2., bearing - 90)[::-1]
        c4 = geo.ll_shift(t_l[1], t_l[0], y_ext / 2., bearing + 90)[::-1]
        c4, c1, c3, c2 = max_width(c4, c1, c3, c2, t_l, t_u)
        # adjust region to fit
        ll_region = (min(c4[0], c1[0], c3[0], c2[0], ll_region[0]), \
                max(c4[0], c1[0], c3[0], c2[0], ll_region[1]), \
                min(c4[1], c1[1], c3[1], c2[1], ll_region[2]), \
                max(c4[1], c1[1], c3[1], c2[1], ll_region[3]))
        with open('/home/vap30/ucgmsim/Velocity-Model/AlpineF2K/Log/VeloModCorners.txt', 'r') as c:
            c.readline()
            c.readline()
            t3 = map(float, c.readline().split())
            t4 = map(float, c.readline().split())
            t1 = map(float, c.readline().split())
            t2 = map(float, c.readline().split())
        print geo.ll_dist(c1[0], c1[1], c4[0], c4[1]), geo.ll_dist(t1[0], t1[1], t4[0], t4[1])
        print geo.ll_dist(c2[0], c2[1], c3[0], c3[1]), geo.ll_dist(t2[0], t2[1], t3[0], t3[1])
        print geo.ll_dist(c1[0], c1[1], c2[0], c2[1]), geo.ll_dist(t1[0], t1[1], t2[0], t2[1])
        print geo.ll_dist(c3[0], c3[1], c4[0], c4[1]), geo.ll_dist(t3[0], t3[1], t4[0], t4[1])

        # proportion in ocean
        ll_region1 = corners2region(c1, c2, c3, c4)
        dxy = auto_dxy(ll_region1)
        path_vm = '%s/landmask_vm.path' % (out)
        grd_vm = '%s/landmask_vm.grd' % (out)
        grd_land = '%s/landmask1.grd' % (out)
        grd_and = '%s/landmask_AND.grd' % (out)
        corners2path(c1, c2, c3, c4, path_vm)
        gmt.grd_mask('f', grd_land, \
                region = ll_region1, dx = dxy, dy = dxy, wd = out)
        gmt.grd_mask(path_vm, grd_vm, \
                region = ll_region1, dx = dxy, dy = dxy, wd = out)
        total_vm = grd_proportion(grd_vm)
        gmt.grdmath([grd_land, grd_vm, 'BITAND', '=', grd_and], \
                region = ll_region1, dx = dxy, dy = dxy, wd = out)
        try:
            land1 = grd_proportion(grd_and) / total_vm * 100
        except ZeroDivisionError:
            # this should not happen, bug in GMT
            # auto_dxy should be a workaround
            land1 = 0

        mid0 = geo.ll_mid(c4[0], c4[1], c2[0], c2[1])
        ylen = geo.ll_dist(c4[0], c4[1], c1[0], c1[1])
        xlen = geo.ll_dist(c4[0], c4[1], c3[0], c3[1])
        if round(land1) < 1:
            xlen = 0
            ylen = 0
    else:
        bearing = 0
        mid0 = geo.ll_mid(min_x[0], min_y[1], max_x[0], max_y[1])
        # following are just x_ext, y_ext??
        ylen = geo.ll_dist(min_x[0], min_x[1], max_x[0], max_x[1])
        xlen = geo.ll_dist(min_y[0], min_y[1], max_y[0], max_y[1])
        if round(land) < 1:
            xlen = 0
            ylen = 0

    # modified sim time
    sim_time_mod = auto_time2(xlen, ylen, 1)

    # assuming hypocentre has been placed at 0.6 * depth
    zmax = auto_z(faultprop.Mw, float(db[dbi + 6].split()[0]) * 0.6)
    # round off
    xlen = round(xlen / hh) * hh
    ylen = round(ylen / hh) * hh
    zlen = round(zmax / hh) * hh

    # store info
    with open('%s/%s.cfg' % (out, name), 'w') as vmd:
        vmd.write('\n'.join(['CALL_TYPE=GENERATE_VELOCITY_MOD', \
                'MODEL_VERSION=1.65', \
                'OUTPUT_DIR=%s' % (name), \
                'ORIGIN_LAT=%s' % (mid0[1]), \
                'ORIGIN_LON=%s' % (mid0[0]), \
                'ORIGIN_ROT=%s' % (bearing), \
                'EXTENT_X=%s' % (xlen), \
                'EXTENT_Y=%s' % (ylen), \
                'EXTENT_ZMAX=%s' % (zlen), \
                'EXTENT_ZMIN=0.0', \
                'EXTENT_Z_SPACING=%s' % (hh), \
                'EXTENT_LATLON_SPACING=%s' % (hh), \
                'MIN_VS=0.5', \
                'TOPO_TYPE=BULLDOZED\n']))
    with open('%s/%s.py' % (out, name), 'w') as pv:
        pv.write('\n'.join(['mag = "%s"' % (faultprop.Mw), \
                'centroidDepth = "%s"' % (float(db[dbi + 6].split()[0]) * 0.6), \
                'MODEL_LAT = "%s"' % (mid0[1]), \
                'MODEL_LON = "%s"' % (mid0[0]), \
                'MODEL_ROT = "%s"' % (bearing), \
                'hh = "%s"' % (hh), \
                'min_vs = "0.5"', \
                'model_version = "1.65"', \
                'topo_type = "BULLDOZED"', \
                'output_directory = "%s"' % (name), \
                'extracted_slice_parameters_directory = "SliceParametersNZ/SliceParametersExtracted.txt"', \
                'code = "rt"', \
                'extent_x = "%s"' % (xlen), \
                'extent_y = "%s"' % (ylen), \
                'extent_zmax = "%s"' % (zlen), \
                'extent_zmin = "0.0"', \
                'sim_duration = "%s"' % (sim_time * (not adjusted) + sim_time_mod * adjusted), \
                'flo = "1.0"', \
                'nx = "%s"' % (int(round(xlen / hh))), \
                'ny = "%s"' % (int(round(ylen / hh))), \
                'nz = "%s"' % (int(round(zlen / hh))), \
                'sufx = "_rt01-h%.3f"' % (hh)]))
    with open(table, 'a') as t:
        if adjusted:
            t.write('%s,%s,%s,%s,%s,%s,%s,%.0f,%s,%s,%s,%s,%.0f\n' % (name, faultprop.Mw, \
                    db[dbi + 6].split()[0], zlen, sim_time, \
                    round(x_ext / hh) * hh, round(y_ext / hh) * hh, \
                    land, zlen, sim_time_mod, xlen, ylen, land1))
        else:
            t.write('%s,%s,%s,%s,%s,%s,%s,%.0f,%s,%s,%s,%s,%.0f\n' % (name, faultprop.Mw, \
                    db[dbi + 6].split()[0], zlen, sim_time, round(x_ext / hh) * hh, round(y_ext / hh) * hh, land, zlen, sim_time_mod, xlen, ylen, land))

    # plot
    p = gmt.GMTPlot('%s/%s.ps' % (out, name))
    p.spacial('M', ll_region, sizing = 7)
    p.coastlines()
    # filled slip area
    p.path('\n'.join([' '.join(map(str, ll)) for ll in pts]), is_file = False, close = True, fill = 'yellow', split = '-')
    # top edge
    p.path('\n'.join([' '.join(map(str, ll)) for ll in pts[:len(pts) / 2]]), is_file = False)
    # vm domain
    #if not adjusted:
    p.path('%s %s\n%s %s\n%s %s\n%s %s' % (min_x[0], min_y[1], max_x[0], min_y[1], max_x[0], max_y[1], min_x[0], max_y[1]), is_file = False, close = True, fill = 'black@95')
    if adjusted:
        p.path('%s %s\n%s %s\n%s %s\n%s %s' % (c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]), is_file = False, close = True, fill = 'black@95', \
                split = '-', width = '1.0p')
    # info text
    p.text((ll_region[0] + ll_region[1]) / 2., ll_region[3], \
            'Mw: %s X: %.0fkm, Y: %.0fkm, land: %.0f%%' \
            % (faultprop.Mw, x_ext, y_ext, land), align = 'CT', dy = -0.1, \
            box_fill = 'white@50')
    if adjusted:
        p.text((ll_region[0] + ll_region[1]) / 2., ll_region[3], \
                'MODIFIED land: %.0f%%' \
                % (land1), align = 'CT', dy = -0.25, \
                box_fill = 'white@50')
    # for testing, show land mask cover
    #p.overlay('90/landmask.grd', 'hot', dx = '2k', dy = '2k', custom_region = ll_region0)
    p.path('res/rough_land.txt', is_file = True, close = False, colour = 'blue', width = '0.2p')
    p.path('res/centre.txt', is_file = True, close = False, colour = 'red', width = '0.2p')
    # actual corners
    p.points('/home/vap30/ucgmsim/Velocity-Model/AlpineF2K/Log/VeloModCorners.txt', fill = 'red', line = None, shape = 'c', size = 0.05)
    p.finalise()
    p.png(dpi = 200, clip = True, background = 'white')

    # clean
    for tmp in ['gmt.conf', 'gmt.history', 'tempEXT.tmp', 'srf.path', \
            'landmask_AND.grd', 'landmask_vm.grd', 'landmask_vm.path', \
            'landmask.grd', 'landmask1.grd', '%s.ps' % (name)]:
        if os.path.exists(os.path.join(out, tmp)):
            os.remove(os.path.join(out, tmp))

    # move to next definition
    dbi += 13 + n_pt

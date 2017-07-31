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

nhm = 'NZ_FLTmodel_2010.txt'
table = 'vminfo.txt'
pgv_target = 2
out = os.path.abspath('out')
centre = 'res/centre.txt'

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
class faultprop:
    Mw = None
    Ztor = 0
    rake = None
    dip = None

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

# proportion of gmt mask file filled (0 -> 1)
def grd_proportion(grd_file):
    try:
        with h5.File(grd_file, 'r') as grd:
            return 100.0 * np.nansum(grd['z'][...]) / float(grd['z'].size)
    except IOError:
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
    over_e = None
    over_w = None
    over_s = None
    for x in np.linspace(0, 1, 81, endpoint = True):
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
        geo.path_from_corners(corners = [ap, bp], output = 'temp000.000', \
                min_edge_points = 80, close = False)
        isections = gmt.intersections(['srf.path', 'res/rough_land.txt', 'temp000.000'], containing = 'temp000.000')
        if len(isections) == 0:
            if over_e == None:
                # land has not been discovered yet
                over_s = x * len_m - 15
            continue
        east = isections[np.argmax(isections, axis = 0)[0]]
        west = isections[np.argmin(isections, axis = 0)[0]]
        #
        m2e = geo.ll_dist(m[0], m[1], east[0], east[1])
        be = geo.ll_bearing(m[0], m[1], east[0], east[1])
        if abs(be - (bearing + 90) % 360) > 90:
            m2e *= -1
        over_e = max(over_e, m2e - m2ee + 15)
        #
        m2w = geo.ll_dist(m[0], m[1], west[0], west[1])
        bw = geo.ll_bearing(m[0], m[1], west[0], west[1])
        if abs(bw - (bearing - 90) % 360) > 90:
            m2w *= -1
        over_w = max(over_w, m2w - m2we + 15)

    # shift
    if over_e < 0:
        a0 = geo.ll_shift(a0[1], a0[0], over_e, bearing + 90)[::-1]
        a1 = geo.ll_shift(a1[1], a1[0], over_e, bearing + 90)[::-1]
    if over_w < 0:
        b0 = geo.ll_shift(b0[1], b0[0], over_w, bearing - 90)[::-1]
        b1 = geo.ll_shift(b1[1], b1[0], over_w, bearing - 90)[::-1]
    if over_s != None:
        a0 = geo.ll_shift(a0[1], a0[0], over_s, bearing_a)[::-1]
        b0 = geo.ll_shift(b0[1], b0[0], over_s, bearing_b)[::-1]
    return a0, a1, b0, b1

# make sure output directory exists for images
if not os.path.exists(out):
    os.makedirs(out)

# clear table file
with open(table, 'w') as _:
    pass

# loop through faults
with open(nhm, 'r') as nf:
    db = nf.readlines()
    dbi = 15
    dbl = len(db)
while dbi < dbl:
    name = db[dbi].strip()
    n_pt = int(db[dbi + 11])
    #if name != 'HikHBaymin':
    #    dbi += 13 + n_pt
    #    continue
    faultprop.Mw = float(db[dbi + 10].split()[0])
    faultprop.rake = float(db[dbi + 5])
    faultprop.dip = float(db[dbi + 3].split()[0])
    rrup, pga_actual = find_rrup()

    # fault width (along dip), projected width (along dip direction)
    fwid = (float(db[dbi + 6].split()[0]) - float(db[dbi + 7].split()[0])) \
            / math.sin(math.radians(faultprop.dip))
    pwid = fwid * math.sin(math.radians(faultprop.dip))
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

    # proportion in ocean
    gmt.grd_mask('f', '%s/landmask.grd' % (out), region = ll_region0, dx = '1k', dy = '1k', wd = out, outside = 0)
    land = grd_proportion('%s/landmask.grd' % (out))

    with open('srf.path', 'w') as sp:
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
        bearing = geo.ll_bearing(mid[0], mid[1], l2, max_y[1])
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

        # proportion in ocean
        ll_region1 = corners2region(c1, c2, c3, c4)
        path_vm = '%s/landmask_vm.path' % (out)
        grd_vm = '%s/landmask_vm.grd' % (out)
        grd_land = '%s/landmask1.grd' % (out)
        grd_and = '%s/landmask_AND.grd' % (out)
        corners2path(c1, c2, c3, c4, path_vm)
        gmt.grd_mask('f', grd_land, \
                region = ll_region1, dx = '1k', dy = '1k', wd = out)
        gmt.grd_mask(path_vm, grd_vm, \
                region = ll_region1, dx = '1k', dy = '1k', wd = out)
        total_vm = grd_proportion(grd_vm)
        gmt.grdmath([grd_land, grd_vm, 'BITAND', '=', grd_and], \
                region = ll_region1, dx = '1k', dy = '1k', wd = out)
        land1 = grd_proportion(grd_and) / total_vm * 100

    # store info
    with open(table, 'a') as t:
        if adjusted:
            t.write('%s,%s,%s,%.0f,%.0f\n' % (name, x_ext, y_ext, land, land1))
        else:
            t.write('%s,%s,%s,%.0f,NaN\n' % (name, x_ext, y_ext, land))

    # plot
    p = gmt.GMTPlot('out/%s.ps' % (name))
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
    #p.overlay('out/west.grd', 'hot', dx = '1k', dy = '1k', custom_region = ll_region0)
    p.path('res/rough_land.txt', is_file = True, close = False, colour = 'blue', width = '0.2p')
    p.path('res/centre.txt', is_file = True, close = False, colour = 'red', width = '0.2p')
    p.finalise()
    p.png(dpi = 200, clip = True, background = 'white')

    # move to next definition
    dbi += 13 + n_pt

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
out = 'out'

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
    #if os.path.exists('%s/%s.png' % (out, name)):
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

    # store info
    with open(table, 'a') as t:
        t.write('%s,%s,%s,%.0f\n' % (name, x_ext, y_ext, land))
    # move to next definition
    dbi += 13 + n_pt
    continue

    # TODO: BELOW

    # modify if necessary
    if faultprop.Mw >= 6.2 and land < 90:
        print 'modifying', name
        gmt.grd_mask('east.txt', '%s/east.grd' % (out), region = ll_region0, dx = '1k', dy = '1k', wd = 'out', outside = 0)
        gmt.grd_mask('west.txt', '%s/west.grd' % (out), region = ll_region0, dx = '1k', dy = '1k', wd = 'out', outside = 0)

    # plot
    p = gmt.GMTPlot('out/%s.ps' % (name))
    p.spacial('M', ll_region, sizing = 7)
    p.coastlines()
    # filled slip area
    p.path('\n'.join([' '.join(map(str, ll)) for ll in pts]), is_file = False, close = True, fill = 'yellow', split = '-')
    # top edge
    p.path('\n'.join([' '.join(map(str, ll)) for ll in pts[:len(pts) / 2]]), is_file = False)
    # vm domain
    p.path('%s %s\n%s %s\n%s %s\n%s %s' % (min_x[0], min_y[1], max_x[0], min_y[1], max_x[0], max_y[1], min_x[0], max_y[1]), is_file = False, close = True, fill = 'black@95')
    # Mw text
    p.text((max_x[0] + min_x[0]) / 2., max_y[1], 'Mw: %s X: %.0fkm, Y: %.0fkm, land: %s' % (faultprop.Mw, x_ext, y_ext, land), align = 'CB', dy = 0.1)
    # for testing, show land mask cover
    p.overlay('out/west.grd', 'hot', dx = '1k', dy = '1k', custom_region = ll_region0)
    p.finalise()
    p.png(dpi = 200, clip = True, background = 'white')

    # move to next definition
    dbi += 13 + n_pt

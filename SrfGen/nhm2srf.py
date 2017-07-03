#!/usr/bin/env python2

from math import sin, radians
import os
import sys

sys.path.append('/home/nesi00213/Pre-processing/SrfGen')
sys.path.append('/home/vap30/ucgmsim/qcore')
sys.path.append('/home/nesi00213/qcore')
from createSRF import CreateSRF_multi
import geo

NHM_FILE = 'NZ_Flt.txt'
GMT_FILE = 'fault_traces.gmt'
N_HYPO = 3

###
### PREPARE
###
if os.path.exists(GMT_FILE):
    os.remove(GMT_FILE)

###
### LOAD FAULTS
###
dbi = 0
with open(NHM_FILE, 'r') as dbr:
    # skip header
    for _ in xrange(15):
        dbr.readline()
    # store rest
    db = map(str.strip, dbr.readlines())
dbl = len(db)

###
### PROCESS FAULTS
###
while dbi < dbl:
    # points in fault trace
    n_pt = int(db[dbi + 11])
    # definition lines to shift to next definition
    n_ln = 13 + n_pt
    pts = [map(float, ll.split()) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
    # clean points (remove duplicates)
    for i in xrange(n_pt - 2, -1, -1):
        if geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) \
                < 0.1:
            del pts[i + 1]
            n_pt -= 1
    # total planes
    n_plane = n_pt - 1
    # strike orientation is based on dip direction
    dip_dir = float(db[dbi + 4])
    strike_avg = (dip_dir - 90) % 360
    # midpoints, lengths and strikes of all faults
    mids = []
    lengths = []
    strikes = []
    for plane in xrange(n_pt - 1):
        mids.append(geo.ll_mid(pts[plane][0], pts[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1]))
        lengths.append(geo.ll_dist(pts[plane][0], pts[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1]))
        bearing = geo.ll_bearing(mids[plane][0], mids[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1])
        if abs(bearing - strike_avg) < 90:
            strikes.append(bearing)
            reverse = False
        else:
            strikes.append((bearing + 180) % 360)
            reverse = True
    # assuming no reverse angles and ordering in one direction
    if reverse:
        mids = mids[::-1]
        lengths = lengths[::-1]
        strikes = strikes[::-1]
    trace_length = sum(lengths)
    hyp_step = trace_length / (N_HYPO * 2.)

    # named values
    seed = 2342
    dt = 0.025
    cases = ['combined']
    name = db[dbi]
    nseg = [n_plane]
    seg_delay = [0]
    mag = [float(db[dbi + 10].split()[0])]
    mom = [-1]
    rvfac_seg = ['-1']
    gwid = ['-1']
    rup_delay = [0]

    rake = [[float(db[dbi + 5])] * n_plane]
    dip = [[float(db[dbi + 3].split()[0])] * n_plane]
    dtop = [[float(db[dbi + 7].split()[0])] * n_plane]
    flen = [lengths]
    dlen = [[0.1] * n_plane]
    fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0]) \
            / sin(radians(dip[0][0]))] * n_plane]
    dwid = [[0.1] * n_plane]
    stk = [strikes]
    elon = [[ll[0] for ll in mids]]
    elat = [[ll[1] for ll in mids]]
    dhypo = [[fwid[0][0] / 2.] * n_plane]

    for n_shyp in xrange(N_HYPO):
        # hypocentre position from far left edge
        shyp_shift = hyp_step * (1 + 2 * n_shyp)
        # NOTE: this shypo is relative to the first combined fault
        # if not adjusted later, must be relative to full length
        shypo = [[shyp_shift - (lengths[0] / 2.)]]
        prefix = '%s_HYP%.2d-%.2d' % (name, n_shyp + 1, N_HYPO)
        # create SRF from description
        CreateSRF_multi(nseg, seg_delay, mag, mom, rvfac_seg, gwid, \
                rup_delay, flen, dlen, fwid, dwid, dtop, stk, rake, dip, \
                elon, elat, shypo, dhypo, dt, seed, prefix, cases)

    # store fault traces
    with open(GMT_FILE, 'a') as traces:
        traces.write('> %s\n%s\n' \
                % (db[dbi], '\n'.join(db[dbi + 12 : dbi + 12 + n_pt])))

    # move to next fault definition
    dbi += n_ln

#!/usr/bin/env python2

from math import sin, radians
import os
import sys

sys.path.append('/home/nesi00213/Pre-processing/SrfGen')
sys.path.append('/home/vap30/ucgmsim/qcore')
sys.path.append('/home/nesi00213/qcore')
from createSRF import CreateSRF_multi
import geo

NHM_FILE = 'NZ_FLTmodel_2010.txt'
GMT_FILE = 'fault_traces.gmt'
N_HYPO = 1

###
### PREPARE
###
if os.path.exists(GMT_FILE):
    os.remove(GMT_FILE)

if len(sys.argv) > 1:
    if not os.path.exists(sys.argv[1]):
        print('Fault selecion file not found: %s' % (sys.argv[1]))
        exit(1)
    with open(sys.argv[1], 'r') as fn_file:
        fault_names = map(str.strip, fn_file.readlines())
else:
    fault_names = None

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
    name = db[dbi]
    # points in fault trace
    n_pt = int(db[dbi + 11])
    # definition lines to shift to next definition
    n_ln = 13 + n_pt
    # skip if not wanted
    if fault_names != None and name not in fault_names:
        dbi += n_ln
        continue
    # load trace
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
    stk_norm = 0
    stk_rev = 0
    for plane in xrange(n_pt - 1):
        mids.append(geo.ll_mid(pts[plane][0], pts[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1]))
        lengths.append(geo.ll_dist(pts[plane][0], pts[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1]))
        bearing = geo.ll_bearing(mids[plane][0], mids[plane][1], \
                pts[plane + 1][0], pts[plane + 1][1])
        if abs((bearing - strike_avg + 180) % 360 - 180) < 90:
            strikes.append(bearing)
            stk_norm += 1
        else:
            strikes.append((bearing + 180) % 360)
            stk_rev += 1
    if stk_norm and stk_rev:
        print('WARNING: FAULT GOES BACK IN REVERSE OF ORIGINAL DIRECTION: %s' \
                % (name))
    # assuming no reverse angles and ordering in one direction
    if stk_rev:
        mids = mids[::-1]
        lengths = lengths[::-1]
        strikes = strikes[::-1]
    trace_length = sum(lengths)
    hyp_step = trace_length / (N_HYPO * 2.)

    # named values
    seed = 2342
    dt = 0.025
    cases = ['combined']
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
    if trace_length < 50:
        dlen = [[0.1] * n_plane]
    else:
        dlen = [[0.2] * n_plane]
    fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0]) \
            / sin(radians(dip[0][0]))] * n_plane]
    dwid = dlen
    stk = [strikes]
    elon = [[ll[0] for ll in mids]]
    elat = [[ll[1] for ll in mids]]
    dhypo = [[fwid[0][0] * 0.6] * n_plane]

    for n_shyp in xrange(N_HYPO):
        # hypocentre position from far left edge
        shyp_shift = hyp_step * (1 + 2 * n_shyp)
        # NOTE: this shypo is relative to the first combined fault
        # if not adjusted later, must be relative to full length
        shypo = [[shyp_shift - (lengths[0] / 2.)]]
        prefix = 'Srf/%s_HYP%.2d-%.2d' % (name, n_shyp + 1, N_HYPO)
        # create SRF from description
        CreateSRF_multi(nseg, seg_delay, mag, mom, rvfac_seg, gwid, \
                rup_delay, flen, dlen, fwid, dwid, dtop, stk, rake, dip, \
                elon, elat, shypo, dhypo, dt, seed, prefix, cases)

    # store fault traces
    #with open(GMT_FILE, 'a') as traces:
    #    traces.write('> %s\n%s\n' \
    #            % (db[dbi], '\n'.join(db[dbi + 12 : dbi + 12 + n_pt])))

    # move to next fault definition
    dbi += n_ln

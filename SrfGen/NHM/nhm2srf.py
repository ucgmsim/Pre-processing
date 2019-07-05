#!/usr/bin/env python2
"""
REQUIREMENTS:
PATH contains plot_srf_square.py plot_srf_map.py if using -p (plot results)
PYTHONPATH contains createSRF and qcore
"""
import sys

import math
import os
from subprocess import call
from time import time

import numpy as np

from createSRF import leonard, skarlatoudis, CreateSRF_multi, gen_meta
from qcore import geo, simulation_structure

FILE_LOCATION = os.path.dirname(os.path.abspath(__file__))

MAGNITUDE_ROUNDING_THRESHOLD = 7.5


def rand_shyp_dhyp(length=1.0, width=1.0):
    # normal distribution
    shyp_mu = 0.5
    shyp_sigma = 0.25
    # weibel distribution
    dhyp_scale = 0.612
    dhyp_shape = 3.353

    shyp = -1
    dhyp = 2

    # within range 0 -> 1 (exclusive)
    while shyp <= 0 or shyp >= 1:
        shyp = np.random.normal(shyp_mu, shyp_sigma)
    while dhyp >= 1:
        dhyp = np.random.weibull(dhyp_shape) * dhyp_scale
    return shyp * length, dhyp * width


def nhm2srf(fault_name, rel_number, out_dir, seed, nhm_skip, nhm_file, plot):
    # adjust outputs
    out_gmt = os.path.join(out_dir, 'fault_traces.gmt')
    out_log = os.path.join(out_dir, 'logfile.txt')
    # prepare file system
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(out_gmt):
        os.remove(out_gmt)
    if not os.path.exists(out_log):
        with open(out_log, 'w') as log:
            log.write('filename\tshyp\tdhyp\tseed\n')

    # load db
    dbi = nhm_skip
    with open(nhm_file, 'r') as dbr:
        db = list(map(str.strip, dbr.readlines()))
    dbl = len(db)

    # Read through the NHM file to find the details for our fault
    while dbi < dbl and db[dbi] != fault_name:
        dbi += 13 + int(db[dbi + 11])

    # If we didn't find it return false
    if db[dbi] != fault_name:
        print("Didn't find fault in nhm")
        return False

    # Load data from NHM
    n_pt = int(db[dbi + 11])
    tect_type, fault_type = db[dbi + 1].split()
    # load trace
    traces = db[dbi + 12: dbi + 12 + n_pt]
    pts = [list(map(float, ll.split())) for ll in traces]
    # clean points (remove duplicates)
    for i in range(n_pt - 2, -1, -1):
        if geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1])  < 0.1:
            del pts[i + 1]
            n_pt -= 1
    # total planes
    n_plane = n_pt - 1

    # other values from nhm file
    mag = [float(db[dbi + 10].split()[0])]
    rake = [[float(db[dbi + 5])] * n_plane]
    dip = [[float(db[dbi + 3].split()[0])] * n_plane]
    dtop = [[float(db[dbi + 7].split()[0])] * n_plane]
    dbottom = float(db[dbi + 6].split()[0])

    # strike orientation is based on dip direction
    dip_dir = float(db[dbi + 4])
    strike_avg = (dip_dir - 90) % 360
    # Finished loading from NHM

    # midpoints, lengths and strikes of all faults
    mids = []
    lengths = []
    strikes = []
    stk_norm = 0
    stk_rev = 0

    for plane in range(n_pt - 1):
        mids.append(geo.ll_mid(pts[plane][0], pts[plane][1], pts[plane + 1][0], pts[plane + 1][1]))
        dist = geo.ll_dist(pts[plane][0], pts[plane][1], pts[plane + 1][0], pts[plane + 1][1])
        lengths.append(round_subfault_size(dist, mag))
        bearing = geo.ll_bearing(mids[plane][0], mids[plane][1], pts[plane + 1][0], pts[plane + 1][1])
        if abs((bearing - strike_avg + 180) % 360 - 180) < 90:
            strikes.append(bearing)
            stk_norm += 1
        else:
            strikes.append((bearing + 180) % 360)
            stk_rev += 1

    # Check direction of the fault
    if stk_norm and stk_rev:
        print('WARNING: FAULT GOES BACK IN REVERSE OF ORIGINAL DIRECTION: {}'.format(fault_name))
        return False

    # assuming no reverse angles and ordering in one direction
    if stk_rev:
        mids = mids[::-1]
        lengths = lengths[::-1]
        strikes = strikes[::-1]
    trace_length = sum(lengths)

    # fixed values
    dt = 0.025
    cases = ['combined']
    nseg = [n_plane]
    seg_delay = [0]
    mom = [-1]
    rvfac_seg = ['-1']
    gwid = ['-1']
    rup_delay = [0]

    flen = [lengths]
    # subfault density
    dlen = [[0.1] * n_plane]
    # Karim: add 3km to depth if bottom >= 12km
    if dbottom >= 12:
        raw_fwid = [[(dbottom - dtop[0][0] + 3) / math.sin(math.radians(dip[0][0]))] * n_plane]
    else:
        raw_fwid = [[(dbottom - dtop[0][0]) / math.sin(math.radians(dip[0][0]))] * n_plane]

    fwid = [[round_subfault_size(f, mag) for f in segment] for segment in raw_fwid]

    # Scale magnitude based on type of interface
    if tect_type == "SUBDUCTION_INTERFACE":
        mag = [skarlatoudis(fwid[0][0] * trace_length)]
    else:
        mag = [leonard(rake[0][0], fwid[0][0] * trace_length)]

    stk = [strikes]
    elon = [[ll[0] for ll in mids]]
    elat = [[ll[1] for ll in mids]]

    shyp_shift, dhyp_shift = rand_shyp_dhyp(trace_length, fwid[0][0])

    shypo = [[shyp_shift - (lengths[0] / 2.)]]
    dhypo = [[dhyp_shift] * n_plane]
    prefix = os.path.join(
        out_dir,
        simulation_structure.get_srf_location(
            simulation_structure.get_realisation_name(
                fault_name, rel_number
            )
        ),
    )[:-4]

    # create SRF from description
    stoch = '{}/{}/Stoch'.format(out_dir, fault_name)
    # store parameters
    with open(out_log, 'a') as log:
        log.write('{}.srf\t{}\t{}\t{}\n'.format(prefix, shypo[0][0], dhypo[0][0], seed))

    # store fault traces
    with open(out_gmt, 'a') as traces:
        traces.write('> {}\n{}\n'.format(fault_name, '\n'.join(traces)))

    t0 = time()
    print('creating SRF: {}'.format(fault_name))
    # all of the work, rest of the script is complete under 1 second
    CreateSRF_multi(nseg, seg_delay, mag, mom, rvfac_seg, gwid, rup_delay, flen, dlen, fwid, dlen, dtop, stk, rake,
                    dip, elon, elat, shypo, dhypo, dt, seed, prefix, cases, dip_dir=dip_dir, stoch=stoch,
                    tect_type=tect_type, silent=True)
    print('created SRF: {} {:.2f}s)'.format(fault_name, time() - t0))

    if plot:
        t0 = time()
        print('plotting SRF: {}'.format(fault_name))
        call(['plot_srf_square.py', '{}.srf'.format(prefix)])
        call(['plot_srf_map.py', '{}.srf'.format(prefix)])
        print('plotted SRF: {} ({:.2f}s)'.format(fault_name, time() - t0))

    srf_file = os.path.join(out_dir, fault_name, 'Srf', "{}_REL01.srf".format(fault_name))

    gen_meta(
        srf_file, 4, mag, stk, rake, dip, dt,
        vm='{}/lp_generic1d-gp01_v1.vmod'.format(os.path.dirname(os.path.abspath(__file__))),
        dip_dir=dip_dir,
        shypo=[s[0] + 0.5 * flen[len(cases) - 1][0] for s in shypo],
        dhypo=[d[0] for d in dhypo], tect_type=tect_type,
        file_name=os.path.join(out_dir, fault_name, fault_name)
    )

    return True


def round_subfault_size(dist, mag):
    if mag[0] > MAGNITUDE_ROUNDING_THRESHOLD:
        return round(dist * 2) / 2
    else:
        return round(dist * 10) / 10


def main():
    from argparse import ArgumentParser

    # parameters
    parser = ArgumentParser()
    parser.add_argument('fault', help='The name of the fault to be generated')
    parser.add_argument(
        'realisation', nargs="?", default=1, type=int, help="The number of the realisation to be generated"
    )
    parser.add_argument('-o', '--out_dir', help='directory to place outputs', default='Sources')
    parser.add_argument(
        '--nhm_file', help='NHM file location', default=os.path.join(FILE_LOCATION, 'NZ_FLTmodel_2010.txt')
    )  # This file is a symbolic link: To avoid need to edit hard-coded filename
    parser.add_argument('--nhm_skip', help='NHM header lines to skip', type=int, default=15)
    parser.add_argument('--dhypo', default=[0.6], type=float, help='depth of hypocentre')
    parser.add_argument('-p', '--plot', help='plot results', action='store_true')
    args = parser.parse_args()
    out_dir = os.path.abspath(args.out_dir)

    fault_name = args.fault

    # load wanted fault information
    if not nhm2srf(fault_name, args.realisation, out_dir, args.seed, args.nhm_skip, args.nhm_file, args.plot):
        sys.exit(1)


if __name__ == '__main__':
    main()

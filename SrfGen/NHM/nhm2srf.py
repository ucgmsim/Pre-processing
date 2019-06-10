#!/usr/bin/env python2
"""
REQUIREMENTS:
PATH contains plot_srf_square.py plot_srf_map.py if using -p (plot results)
PYTHONPATH contains createSRF and qcore
"""

import math
import os
from subprocess import call
import sys
from time import time

import numpy as np

from createSRF import leonard, skarlatoudis, CreateSRF_multi
from qcore import geo, simulation_structure


def rand_shyp_dhyp(length = 1.0, width = 1.0):
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

###
### PREPARE TASKS
###
# most work is trivial, only need to run CreateSRF_multi parallel
def load_msgs(args, fault_names, faults):
    # adjust outputs
    out_gmt = os.path.join(args.out_dir, 'fault_traces.gmt')
    out_log = os.path.join(args.out_dir, 'logfile.txt')
    # prepare file system
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
    if os.path.exists(out_gmt):
        os.remove(out_gmt)
    with open(out_log, 'w') as log:
        log.write('filename\tshyp\tdhyp\tseed\n')


    # load db
    dbi = args.nhm_skip
    with open(args.nhm_file, 'r') as dbr:
        db = list(map(str.strip, dbr.readlines()))
    dbl = len(db)

    # compile messages (parameters for CreateSRF_multi)
    msgs = []
    while dbi < dbl:
        name = db[dbi]
        tect_type, fault_type = db[dbi + 1].split()
        # points in fault trace
        n_pt = int(db[dbi + 11])
        skip = 13 + n_pt
        # skip if not wanted
        if fault_names != None and name not in fault_names:
            dbi += skip
            continue

        # load trace
        pts = [list(map(float, ll.split())) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
        # clean points (remove duplicates)
        for i in range(n_pt - 2, -1, -1):
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
        for plane in range(n_pt - 1):
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
            print('WARNING: FAULT GOES BACK '
                    'IN REVERSE OF ORIGINAL DIRECTION: %s' \
                    % (name))
        # assuming no reverse angles and ordering in one direction
        if stk_rev:
            mids = mids[::-1]
            lengths = lengths[::-1]
            strikes = strikes[::-1]
        trace_length = sum(lengths)

        # wanted parameters to override
        t_hypo = 'n'
        n_hypo = args.nhypo
        n_slip = args.nslip
        dhypos = args.dhypo
        if faults != None:
            fault = faults[fault_names.index(name)]
            if len(fault) >= 2:
                # given as hypocentre every x km
                if fault[1][-1] == 'k':
                    t_hypo = 'k'
                    hyp_step = float(fault[1][:-1])
                    n_hypo = 1 + int(trace_length // hyp_step)
                    # 0th hypocentre position
                    if n_hypo == 1:
                        z_hypo = trace_length / 2.
                    else:
                        z_hypo = (trace_length % hyp_step) / 2.

                # given as number of randomly placed hypocentres
                elif fault[1][-1] == 'r':
                    t_hypo = 'r'
                    n_hypo = int(fault[1][:-1])

                # given as number of hypocentres
                else:
                    n_hypo = int(fault[1])

            if len(fault) >= 3:
                n_slip = int(fault[2])
            if len(fault) >= 4:
                dhypos = list(map(float, fault[3].split(',')))
        if t_hypo == 'n':
            hyp_step = trace_length / (n_hypo * 2.)
        seed = args.seed

        # fixed values
        dt = 0.025
        cases = ['combined']
        nseg = [n_plane]
        seg_delay = [0]
        mom = [-1]
        rvfac_seg = ['-1']
        gwid = ['-1']
        rup_delay = [0]

        # other values from nhm file
        mag = [float(db[dbi + 10].split()[0])]
        rake = [[float(db[dbi + 5])] * n_plane]
        dip = [[float(db[dbi + 3].split()[0])] * n_plane]
        dtop = [[float(db[dbi + 7].split()[0])] * n_plane]
        flen = [lengths]
        # subfault density
        dlen = [[0.1] * n_plane]
        fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0]) \
                / math.sin(math.radians(dip[0][0]))] * n_plane]
        # Karim: add 3km to depth if bottom >= 12km
        if float(db[dbi + 6].split()[0]) >= 12:
            fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0] + 3) \
                / math.sin(math.radians(dip[0][0]))] * n_plane]
        if tect_type == "SUBDUCTION_INTERFACE":
            mag = [skarlatoudis(fwid[0][0] * trace_length)]
        else:
            mag = [leonard(rake[0][0], fwid[0][0] * trace_length)]
        dwid = dlen
        stk = [strikes]
        elon = [[ll[0] for ll in mids]]
        elat = [[ll[1] for ll in mids]]

        for n_shyp in range(n_hypo):
            # hypocentre position from far left edge
            if t_hypo == 'n':
                shyp_shift = hyp_step * (1 + 2 * n_shyp)
            elif t_hypo == 'k':
                shyp_shift = z_hypo + hyp_step * n_shyp
            elif t_hypo == 'r':
                shyp_shift, dhyp_shift = rand_shyp_dhyp(trace_length, fwid[0][0])
            # NOTE: this shypo is relative to the first combined fault
            # if not adjusted later, must be relative to full length
            shypo = [[shyp_shift - (lengths[0] / 2.)]]
            for _ in range(n_slip):
                for i, d in enumerate(dhypos):
                    seed += args.seed_inc
                    if t_hypo == 'r':
                        dhypo = [[dhyp_shift] * n_plane]
                    else:
                        dhypo = [[fwid[0][0] * d] * n_plane]
                    prefix = os.path.join(
                        args.out_dir,
                        simulation_structure.get_srf_location(
                            simulation_structure.get_realisation_name(
                                name, n_shyp+1
                            )
                        ),
                    )[:-4]
                    # create SRF from description
                    msgs.append({'nseg':nseg, 'seg_delay':seg_delay, 'mag':mag, \
                            'mom':mom, 'rvfac_seg':rvfac_seg, 'gwid':gwid, \
                            'rup_delay':rup_delay, 'flen':flen, 'dlen':dlen, \
                            'fwid':fwid, 'dwid':dwid, 'dtop':dtop, 'stk':stk, \
                            'rake':rake, 'dip':dip, 'elon':elon, 'elat':elat, \
                            'shypo':shypo, 'dhypo':dhypo, 'dt':dt, 'seed':seed, \
                            'prefix':prefix, 'cases':cases, 'dip_dir':dip_dir, \
                            'stoch':'%s/%s/Stoch' % (args.out_dir, name), \
                            'name':name, 'tect_type':db[dbi + 1].split()[0], \
                            'plot':args.plot})
                    # store parameters
                    with open(out_log, 'a') as log:
                        log.write('%s.srf\t%s\t%s\t%s\n' \
                                % (prefix, shypo[0][0], dhypo[0][0], seed))

        # store fault traces
        with open(out_gmt, 'a') as traces:
            traces.write('> %s\n%s\n' \
                    % (name, '\n'.join(db[dbi + 12 : dbi + 12 + n_pt])))

        # move to next fault definition
        dbi += skip

    return msgs

def run_create_srf(t):
    t0 = time()
#    sys.stdout = open(str(os.getpid())+".out","w")
    print('creating SRF: %s' % (t['name']))
    # all of the work, rest of the script is complete under 1 second
    CreateSRF_multi(t['nseg'], t['seg_delay'], t['mag'], t['mom'], \
            t['rvfac_seg'], t['gwid'], t['rup_delay'], t['flen'], t['dlen'], \
            t['fwid'], t['dwid'], t['dtop'], t['stk'], t['rake'], t['dip'], \
            t['elon'], t['elat'], t['shypo'], t['dhypo'], t['dt'], t['seed'], \
            t['prefix'], t['cases'], dip_dir = t['dip_dir'], \
            stoch = t['stoch'], tect_type = t['tect_type'], silent = True)
    print('created SRF: %s (%.2fs)' % (t['name'], time() - t0))
    if t['plot']:
        t0 = time()
        print('plotting SRF: %s' % (t['name']))
        call(['plot_srf_square.py', '%s.srf' % (t['prefix'])])
        call(['plot_srf_map.py', '%s.srf' % (t['prefix'])])
        print('plotted SRF: %s (%.2fs)' % (t['name'], time() - t0))
 #   sys.stdout.close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    from multiprocessing import Pool

    # parameters
    parser = ArgumentParser()
    arg = parser.add_argument
    arg('selection_file', help = 'fault selection file')
    parser.add_argument('-o', '--out-dir', help = 'directory to place outputs', \
            default = 'autosrf')
    arg('--nhm-file', help = 'NHM file location', \
        default = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
                               'NZ_FLTmodel_2010.txt')) # This file is a symbolic link: To avoid need to edit hard-coded filename
    arg('--nhm-skip', help = 'NHM header lines to skip', \
            type = int, default = 15)
    arg('--seed', help = 'initial seed', type = int, default = 1234)
    arg('--seed-inc', help = 'seed increment', type = int, default = 10)
    arg('--nhypo', default = 1, \
        help = 'hypocentre number/spacing if not in selection_file')
    arg('--nslip', type = int, default = 1, \
        help = 'number of slips if not in selection_file')
    arg('--dhypo', action = 'append', type = float, \
        help = 'depth of hypocentre, repeat parameter as needed')
    arg('-n', '--nproc', help = 'number of processes', type = int, default = 1)
    arg('-p', '--plot', help = 'plot results', action = 'store_true')
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)
    # selection file or ALL
    if args.selection_file == 'ALL':
        faults = None
        fault_names = None
    elif not os.path.exists(args.selection_file):
        print('Fault selecion file not found: %s' % (args.selection_file))
        sys.exit(1)
    else:
        with open(args.selection_file, 'r') as select:
            faults = list(map(str.split, select.readlines()))
        fault_names = [f[0] for f in faults]
    # default value
    if args.dhypo is None:
        args.dhypo = [0.6]

    # load wanted fault information
    msg_list = load_msgs(args, fault_names, faults)
    if len(msg_list) == 0:
        print('No matches found.')
        sys.exit(1)

    # distribute work
    p = Pool(args.nproc)
    p.map(run_create_srf, msg_list)
    # debug friendly alternative
    #[run_create_srf(msg) for msg in msg_list]

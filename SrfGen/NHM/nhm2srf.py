#!/usr/bin/env python2

from math import sin, radians, log10
import os
import sys

from mpi4py import MPI

# .. points to createSRF
sys.path.append('..')
from createSRF import CreateSRF_multi
from qcore import geo

NHM_FILE = 'NZ_FLTmodel_2010.txt'
NHM_START = 15
GMT_FILE = 'fault_traces.gmt'
LOG_FILE = 'logfile.txt'
N_HYPO = 1
N_SLIP = 1
SEED_0 = 1234
SEED_INC = 10
# MPI - do not change
MASTER = 0

# Leonard 2014 Relations
def leonard(rake, A):
    # if dip slip else strike slip
    if round(rake % 360 / 90.) % 2:
        return 4.19 + log10(A)
    else:
        return 4.18 + log10(A)

###
### PREPARE TASKS
###
# most work is trivial, only need to run CreateSRF_multi parallel
def load_msgs(fault_names, faults, out):
    # adjust outputs
    out_gmt = os.path.join(out, GMT_FILE)
    out_log = os.path.join(out, LOG_FILE)
    # prepare file system
    if not os.path.isdir(out):
        os.makedirs(out)
    if os.path.exists(out_gmt):
        os.remove(out_gmt)
    with open(out_log, 'w') as log:
        log.write('filename\tshyp\tdhyp\tseed\n')


    # load db
    dbi = NHM_START
    with open(NHM_FILE, 'r') as dbr:
        db = map(str.strip, dbr.readlines())
    dbl = len(db)

    # compile messages (parameters for CreateSRF_multi)
    msgs = []
    while dbi < dbl:
        name = db[dbi]
        # points in fault trace
        n_pt = int(db[dbi + 11])
        skip = 13 + n_pt
        # skip if not wanted
        if fault_names != None and name not in fault_names:
            dbi += skip
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
        if fault_names == None:
            n_hypo = N_HYPO
            n_slip = N_SLIP
        else:
            fault = faults[fault_names.index(name)]
            try:
                # given as number of hypocentres
                n_hypo = int(fault[1])
            except ValueError:
                # given as hypocentre every x km
                t_hypo = 'k'
                hyp_step = float(fault[1][:-1])
                n_hypo = 1 + int(trace_length // hyp_step)
                # 0th hypocentre position
                if n_hypo == 1:
                    z_hypo = trace_length / 2.
                else:
                    z_hypo = (trace_length % hyp_step) / 2.
            n_slip = int(fault[2])
        if t_hypo == 'n':
            hyp_step = trace_length / (n_hypo * 2.)
        seed = SEED_0

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
        if trace_length < 50:
            dlen = [[0.1] * n_plane]
        elif trace_length < 100:
            dlen = [[0.2] * n_plane]
        else:
            dlen = [[0.5] * n_plane]
        fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0]) \
                / sin(radians(dip[0][0]))] * n_plane]
        # Karim: add 3km to depth if bottom >= 12km
        if float(db[dbi + 6].split()[0]) >= 12:
            fwid = [[(float(db[dbi + 6].split()[0]) - dtop[0][0] + 3) \
                / sin(radians(dip[0][0]))] * n_plane]
            mag = [leonard(rake[0][0], fwid[0][0] * trace_length)]
        dwid = dlen
        stk = [strikes]
        elon = [[ll[0] for ll in mids]]
        elat = [[ll[1] for ll in mids]]
        dhypo = [[fwid[0][0] * 0.6] * n_plane]

        for n_shyp in xrange(n_hypo):
            # hypocentre position from far left edge
            if t_hypo == 'n':
                shyp_shift = hyp_step * (1 + 2 * n_shyp)
            elif t_hypo == 'k':
                shyp_shift = z_hypo + hyp_step * n_shyp
            # NOTE: this shypo is relative to the first combined fault
            # if not adjusted later, must be relative to full length
            shypo = [[shyp_shift - (lengths[0] / 2.)]]
            for _ in xrange(n_slip):
                seed += SEED_INC
                prefix = '%s/%s/Srf/%s_HYP%.2d-%.2d_S%s' \
                        % (out, name, name, n_shyp + 1, n_hypo, seed)
                # create SRF from description
                msgs.append({'nseg':nseg, 'seg_delay':seg_delay, 'mag':mag, \
                        'mom':mom, 'rvfac_seg':rvfac_seg, 'gwid':gwid, \
                        'rup_delay':rup_delay, 'flen':flen, 'dlen':dlen, \
                        'fwid':fwid, 'dwid':dwid, 'dtop':dtop, 'stk':stk, \
                        'rake':rake, 'dip':dip, 'elon':elon, 'elat':elat, \
                        'shypo':shypo, 'dhypo':dhypo, 'dt':dt, 'seed':seed, \
                        'prefix':prefix, 'cases':cases, 'dip_dir':dip_dir, \
                        'stoch':'%s/%s/Stoch' % (out, name), 'name':name})
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
    # all of the work, rest of the script is complete under 1 second
    CreateSRF_multi(t['nseg'], t['seg_delay'], t['mag'], t['mom'], \
            t['rvfac_seg'], t['gwid'], t['rup_delay'], t['flen'], t['dlen'], \
            t['fwid'], t['dwid'], t['dtop'], t['stk'], t['rake'], t['dip'], \
            t['elon'], t['elat'], t['shypo'], t['dhypo'], t['dt'], t['seed'], \
            t['prefix'], t['cases'], dip_dir = t['dip_dir'], \
            stoch = t['stoch'], silent = True)

###
### MASTER
###
if len(sys.argv) > 1:
    # clock start
    t_start = MPI.Wtime()

    # PARAMETER 1: name of nhm selection file
    if sys.argv[1] == 'ALL':
        faults = None
        fault_names = None
    elif not os.path.exists(sys.argv[1]):
        print('Fault selecion file not found: %s' % (sys.argv[1]))
        sys.exit(1)
    else:
        with open(sys.argv[1], 'r') as select:
            faults = map(str.split, select.readlines())
        fault_names = [f[0] for f in faults]
    # PARAMETER 2: output folder (optional)
    if len(sys.argv) > 2:
        out = os.path.abspath(sys.argv[2])
    else:
        out = os.path.abspath('autosrf')

    # load wanted fault information
    msg_list = load_msgs(fault_names, faults, out)
    if len(msg_list) == 0:
        print('No matches found.')
        sys.exit(1)

    # do not fully load a machine with many cores
    nproc_max = int(os.sysconf('SC_NPROCESSORS_ONLN'))
    if nproc_max > 8:
        nproc_max = int(round(0.8 * nproc_max))
    # not more processes than jobs
    nproc = min(len(msg_list), nproc_max)
    # spawn slaves
    comm = MPI.COMM_WORLD.Spawn(
        sys.executable, args = [sys.argv[0]], maxprocs = nproc)
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

    # gather, reports aren't used
    reports = comm.gather(None, root = MPI.ROOT)
    # stop mpi
    comm.Disconnect()

###
### SLAVE
###
else:
    # connect to parent
    try:
        comm = MPI.Comm.Get_parent()
        rank = comm.Get_rank()
    except MPI.Exception:
        print('First parameter is fault selection file.')
        print('Use "ALL" to loop through all faults.')
        print('First parameter not found.')
        print('Alternatively MPI cannot connect to parent.')
        sys.exit(1)

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):
        # task timer
        t0 = MPI.Wtime()
        # run
        try:
            details = run_create_srf(task)
        except:
            print('FAILED TASK: %s' % (task['name']))
            continue
        # store
        logbook.append(details)

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()

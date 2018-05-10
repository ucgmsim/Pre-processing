#!/usr/bin/env python2

from argparse import ArgumentParser
import os
from subprocess import call
import sys

from mpi4py import MPI
import numpy as np

from createSRF import leonard, CreateSRF_ps

# MPI
MASTER = 0

def run_create_srf(args, t, vs, rho):
    mom = -1
    prefix = os.path.join(args.out_dir, str(t.pid))
    CreateSRF_ps(t.lat, t.lon, t.depth, t.mag, \
            mom, t.strike, t.rake, t.dip, dt = args.dt, \
            prefix = prefix, stoch = args.out_dir, \
            vs = vs, rho = rho, silent = True)
    if args.plot:
        call(['plot_srf_map.py', '%s.srf' % (prefix)])

###
### MASTER
###
if len(sys.argv) > 1:
    # load parameters
    parser = ArgumentParser()
    parser.add_argument('csv_file', help = 'path to CMT solutions CSV')
    parser.add_argument('velocity_model', help = 'path to 1D velocity model')
    parser.add_argument('-o', '--out-dir', help = 'directory to place outputs', \
            default = 'autosrf')
    parser.add_argument('--dt', help = 'SRF timestep', type = float, \
            default = 0.005)
    parser.add_argument('-n', '--nproc', type = int, \
            default = int(os.sysconf('SC_NPROCESSORS_ONLN') - 1))
    parser.add_argument('-p', '--plot', action = 'store_true')
    args = parser.parse_args()
    # validate parameters
    if not os.path.exists(args.csv_file):
        print('CMT solutions CSV file not found: %s' % (args.csv_file))
        sys.exit(1)
    if not os.path.exists(args.velocity_model):
        print('Velocity model file not found: %s' % (args.velocity_model))
        sys.exit(1)

    # load CMT CSV
    sources = np.rec.array(np.loadtxt(args.csv_file, dtype = [('pid', 'i4'), \
            ('lat', 'f4'), ('lon', 'f4'), ('depth', 'f4'), ('mag', 'f4'), \
            ('strike', 'f4'), ('dip', 'f4'), ('rake', 'f4')], \
            delimiter = ',', skiprows = 1, usecols = (0, 2, 3, 13, 11, 4, 5, 6)))

    # load velocity model (vs and rho)
    vmodel = np.rec.array(np.loadtxt(args.velocity_model, \
            skiprows = 1, usecols = (0, 2, 3), \
            dtype = [('depth', 'f8'), ('vs', 'f8'), ('rho', 'f8')]))
    vmodel.depth = np.cumsum(vmodel.depth)
    vs = np.interp(sources.depth, vmodel.depth, vmodel.vs)
    rho = np.interp(sources.depth, vmodel.depth, vmodel.rho)

    # not more processes than jobs
    args.nproc = min(len(sources), args.nproc)
    # spawn slaves
    comm = MPI.COMM_WORLD.Spawn(
        sys.executable, args = [sys.argv[0]], maxprocs = args.nproc)
    status = MPI.Status()

    # distribute work to slaves who ask
    for i, s in enumerate(sources):
        # slave asking for work
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        slave_id = status.Get_source()
        comm.send(obj = [args, s, vs[i], rho[i]], dest = slave_id)
    # kill slaves
    for _ in xrange(args.nproc):
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        slave_id = status.Get_source()
        # kill signal
        comm.send(obj = StopIteration, dest = slave_id)

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
        print('First parameter is CMT solution CSV.')
        print('First parameter not found.')
        print('Alternatively MPI cannot connect to parent.')
        sys.exit(1)

    # ask for work until stop sentinel
    logbook = []
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):
        # task timer
        t0 = MPI.Wtime()
        # run
        details = run_create_srf(task[0], task[1], task[2], task[3])
        # store
        #logbook.append(details)

    # reports to master
    comm.gather(sendobj = logbook, root = MASTER)
    # shutdown
    comm.Disconnect()

#!/usr/bin/env python2

from argparse import ArgumentParser
import os
from subprocess import call
import sys

from mpi4py import MPI
import numpy as np

from createSRF import leonard, CreateSRF_ps

# MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
MASTER = 0
if size < 2:
    sys.exit('use mpirun with nproc > 1')

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
if rank == MASTER:
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
    sources = np.rec.array(np.loadtxt(args.csv_file, dtype = [('pid', '|S32'), \
            ('lat', 'f8'), ('lon', 'f8'), ('depth', 'f8'), ('mag', 'f8'), \
            ('strike', 'f8'), ('dip', 'f8'), ('rake', 'f8')], \
            delimiter = ',', skiprows = 1, usecols = (0, 2, 3, 13, 11, 4, 5, 6)))

    # load velocity model (vs and rho)
    vmodel = np.rec.array(np.loadtxt(args.velocity_model, \
            skiprows = 1, usecols = (0, 2, 3), \
            dtype = [('depth', 'f8'), ('vs', 'f8'), ('rho', 'f8')]))
    depth_bins = np.digitize(np.around(sources.depth, decimals = 5), \
                             np.around(np.cumsum(vmodel.depth), decimals = 5))
    vs = vmodel.vs[depth_bins]
    rho = vmodel.rho[depth_bins]

    # distribute work to slaves who ask
    status = MPI.Status()
    for i, s in enumerate(sources):
        # slave asking for work
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        slave_id = status.Get_source()
        comm.send(obj = [args, s, vs[i], rho[i]], dest = slave_id)
    # kill slaves
    for _ in xrange(size - 1):
        comm.recv(source = MPI.ANY_SOURCE, status = status)
        slave_id = status.Get_source()
        # kill signal
        comm.send(obj = StopIteration, dest = slave_id)

###
### SLAVE
###
else:
    # ask for work until stop sentinel
    for task in iter(lambda: comm.sendrecv(None, dest = MASTER), StopIteration):
        run_create_srf(task[0], task[1], task[2], task[3])

#!/usr/bin/env python2

from argparse import ArgumentParser
from multiprocessing import Pool
import os
from subprocess import call
import sys

import numpy as np

from createSRF import leonard, CreateSRF_ps

def run_create_srf(args, t, vs, rho):
    mom = -1
    prefix = os.path.join(args.out_dir, str(t.pid))
    CreateSRF_ps(t.lat, t.lon, t.depth, t.mag, \
            mom, t.strike, t.rake, t.dip, dt = args.dt, \
            prefix = prefix, stoch = args.out_dir, \
            vs = vs, rho = rho, silent = True)
    if args.plot:
        call(['plot_srf_map.py', '%s.srf' % (prefix)])

def run_create_srf_star(args_t_vs_rho):
    return run_create_srf(*args_t_vs_rho)

# load parameters
parser = ArgumentParser()
parser.add_argument('csv_file', help = 'path to CMT solutions CSV')
parser.add_argument('velocity_model', help = 'path to 1D velocity model')
parser.add_argument('-o', '--out-dir', help = 'directory to place outputs', \
                    default = 'autosrf')
parser.add_argument('--dt', help = 'SRF timestep', type = float, \
                    default = 0.005)
parser.add_argument('-n', '--nproc', type = int, default = 1)
parser.add_argument('-p', '--plot', action = 'store_true')
args = parser.parse_args()
# validate parameters
if not os.path.exists(args.csv_file):
    sys.exit('CMT solutions CSV file not found: %s' % (args.csv_file))
if not os.path.exists(args.velocity_model):
    sys.exit('Velocity model file not found: %s' % (args.velocity_model))

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

# distribute work
p = Pool(args.nproc)
p.map(run_create_srf_star, zip([args] * len(sources), sources, vs, rho))
#[run_create_srf(args, sources[i], vs[i], rho[i]) for i in xrange(len(sources))]

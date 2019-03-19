#!/usr/bin/env python2

from argparse import ArgumentParser
from multiprocessing import Pool
import os
from subprocess import call
import sys

import numpy as np
import pandas as pd
import yaml

from createSRF import leonard, CreateSRF_ps
from createSourceRealisation import create_ps_realisation


def run_create_srf(args, t, vs, rho, n_sims):
    mom = -1
    try:
        pid=t.pid.decode() #t.pid is originally numpy.bytes_
    except AttributeError:
        pid = t.pid

    prefix = os.path.join(args.out_dir, pid, 'Srf', pid)
    if args.uncertainty_file:
        with open(args.uncertainty_file) as ao_f:
            add_opts = yaml.load(ao_f)
        create_ps_realisation(args.out_dir, pid, t.lat, t.lon, t.depth, t.mag, mom, t.strike, t.rake, t.dip,
                              n_realisations=n_sims, additional_options=add_opts, dt=args.dt, vs=vs, rho=rho,
                              silent=True)
    else:
        stoch = os.path.join(args.out_dir, pid, 'Stoch')
        CreateSRF_ps(t.lat, t.lon, t.depth, t.mag, mom, t.strike, t.rake, t.dip, dt=args.dt, prefix=prefix, stoch=stoch,
                     vs=vs, rho=rho, silent=True)
    if args.plot:
        call(['plot_srf_map.py', '%s.srf' % (prefix)])


def run_create_srf_star(args_t_vs_rho):
    return run_create_srf(*args_t_vs_rho)


# load parameters
parser = ArgumentParser()
parser.add_argument('csv_file', help = 'path to CMT solutions CSV')
parser.add_argument('velocity_model', help = 'path to 1D velocity model')
parser.add_argument('-o', '--out-dir', help = 'directory to place outputs', \
                    default = './autosrf')
parser.add_argument('--dt', help = 'SRF timestep', type = float, \
                    default = 0.005)
parser.add_argument('-n', '--nproc', type = int, default = 1)
parser.add_argument('-p', '--plot', action = 'store_true')
uncertainty_group = parser.add_argument_group('Uncertainty options', 'Options to be used when SRFs are to be genreated '
                                                                     'with uncertainty factors')
uncertainty_group.add_argument('-r', '--realistion-uncertainty', dest='uncertainty_file', type=str, default='')
uncertainty_group.add_argument('-c', '--cybershake-fault-file', dest='cs_file', type=str, default='')
args = parser.parse_args()
# validate parameters
if not os.path.exists(args.csv_file):
    sys.exit('CMT solutions CSV file not found: %s' % (args.csv_file))
if not os.path.exists(args.velocity_model):
    sys.exit('Velocity model file not found: %s' % (args.velocity_model))
all_opts = bool(args.uncertainty_file) and bool(args.cs_file)
any_opts = bool(args.uncertainty_file) or bool(args.cs_file) or bool(args.add_opts_file)
if any_opts and not all_opts:
    parser.error("If realisation uncertainty, a cybershake fault file or an additional options file are given then the "
                 "others must be given also.")
    exit(1)


# load CMT CSV
sources = np.rec.array(np.loadtxt(args.csv_file, dtype = [('pid', '|S32'), \
        ('lat', 'f8'), ('lon', 'f8'), ('depth', 'f8'), ('mag', 'f8'), \
        ('strike', 'f8'), ('dip', 'f8'), ('rake', 'f8')], \
        delimiter = ',', skiprows = 1, usecols = (0, 2, 3, 13, 11, 4, 5, 6)))

if all_opts:
    # Sort the array by pid to match the pandas array
    sources = np.sort(sources, order='pid')

    with open(args.cs_file) as cs_f:
        cs_file_pd = pd.read_csv(cs_f, sep=' ', header=None, names=['pid', 'N_Rels'])

    # Apply the identity function to the PublicID column and a filter function to the second column
    cs_file_pd = cs_file_pd.transform({'pid': lambda x: x,
                                       'N_Rels': lambda x: int(x[:-1]) if x.endswith('r') else 0})

    # Remove all values that failed the filter function
    cs_file_pd = cs_file_pd[cs_file_pd.N_Rels != 0]

    n_sims = list(cs_file_pd.sort_values('pid')['N_Rels'])
else:
    n_sims = np.ones(len(sources))

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
p.map(run_create_srf_star,list(zip([args] * len(sources), sources, vs, rho, n_sims)))
#[run_create_srf(args, sources[i], vs[i], rho[i]) for i in xrange(len(sources))]

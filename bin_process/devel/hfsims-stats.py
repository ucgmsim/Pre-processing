#!/usr/bin/env python2
"""
Simulates high frequency seismograms for stations.

USAGE: edit params.py as needed. Then, execute this file.
For parameter documentation, see params.py or the simulator source code.

Original file from Rob (Jan 2015).
Converted to python.
@date 05 May 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
"""
import sys
from multiprocessing import Pool
import os
sys.path.append(os.path.abspath(os.path.curdir))
from subprocess import call, Popen, PIPE

from shared import *
from params import *
from params_base_bb import *

local_statfile = os.path.join(hf_sim_dir, 'local.statfile')
procs = 64
if len(sys.argv) > 1:
    procs = int(sys.argv[1])
if len(sys.argv) > 2:
    hf_sim_bin = os.path.join(global_root, sys.argv[2])

# verify input incl. params.py
verify_binaries([hf_sim_bin])
verify_files([stat_file, hf_slip, hf_v_model])
verify_strings([hf_prefix, hf_t_len, hf_dt, hf_vs_moho, hf_fa_sig_1, hf_rv_sig_1, \
        hf_sdrop, hf_kappa, hf_qfexp, hf_rayset, hf_rvfac, hf_shal_rvfac, hf_deep_rvfac, \
        hf_czero, hf_site_amp, hf_mom, hf_rupv, hf_seed, hf_fmax, hf_calpha])
verify_user_dirs([hf_sim_dir])
verify_logfiles([local_statfile])
verify_user_dirs([hf_accdir],reset=True)

# don't give a single process too many stations to work on
with open(stat_file, 'r') as fp:
    stat_estimate = len(fp.readlines())
if stat_estimate > procs * 500:
    jobs = stat_estimate / 500
else:
    jobs = procs

statfile_base = os.path.join(hf_sim_dir, 'local.tmp.sf_')
statfiles = ['%s%d' % (statfile_base, p) for p in xrange(jobs)]

out_prefix = os.path.join(hf_accdir, hf_prefix)

# copy line not starting with '#' into a local stat_file copy
# file pointers for each processes' stat_file portion
fps = []
# number of stations for each process
nss = [0] * jobs
for f in statfiles:
    fps.append(open(f, 'w'))
with open(stat_file, 'r') as fp:
    # in fp.readline() instead of in fp for AIX/old python? compatibility
    for line in fp.readlines():
        if line[0] != '#':
            fps[sum(nss) % jobs].write(line)
            nss[sum(nss) % jobs] += 1
for fp in fps:
    # file must end with new line
    fp.write('\n')
    fp.close()
# some info echoed in original script, might be useful for some
print("%d %s" % (sum(nss), hf_v_model))

def run_hf((local_statfile, n_stat)):
    # just in case more processes than stations
    if n_stat == 0:
        return
    # prevent terminal echo slow down, not useful with many processes anyway
    sink = open('/dev/null', 'w')
    # execute calculation, stdin from pipe
    hf_pipe = Popen([hf_sim_bin], stdin = PIPE, stdout = sink)
    # pipe input data into binary stdin
    # look at hf_sim_bin source code (same dir) for more info on parameters
    hf_pipe.communicate('\n%s\n%s\n%s\n%s\n%s\n4   0  0.02  19.9\n%s\n%s\n%s %s %s %s %s\n%s %s %s %s %s\n%s %s\n%s\n%s\n%s\n-99 0.0 0.0 0.0 0.0 1\n-1\n%s 0.0 %s\n1\n' % (hf_sdrop, local_statfile, out_prefix, hf_rayset, hf_site_amp, hf_seed, n_stat, hf_t_len, hf_dt, hf_fmax, hf_kappa, hf_qfexp, hf_rvfac, hf_shal_rvfac, hf_deep_rvfac, hf_czero, hf_calpha, hf_mom, hf_rupv, hf_slip, hf_v_model, hf_vs_moho, hf_fa_sig_1, hf_rv_sig_1))
    sink.close()

p = Pool(procs)
p.map(run_hf, zip(statfiles, nss))

for temp_statfile in statfiles:
    os.remove(temp_statfile)

set_permission(hf_sim_dir)

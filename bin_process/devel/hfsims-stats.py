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

#hf_sim_dir = '.'
#hf_sim_bin = './hb_high_mod'
#stat_file = 'cantstations_0.5_h0.100.ll'
#hf_slip = 'm6.20-16.0x9.0_s560.stoch'
#hf_v_model = 'Cant1D_v2-midQ.1d'
#hf_prefix = 'test'
#hf_t_len = '100' # seconds
#hf_dt = '0.005'
#hf_vs_moho = '999.9'
#hf_fa_sig_1 = '0.0'
#hf_rv_sig_1 = '0.1'
#hf_sdrop = '50'
#hf_kappa = '0.045'
#hf_qfexp = '0.6'
#hf_rayset = '2 1 2'
#hf_rvfac = '0.8'
#hf_shal_rvfac = '0.70'
#hf_deep_rvfac = '0.70'
#hf_czero = '2.1'
#hf_site_amp = '1'
#hf_mom = '-1'
#hf_rupv = '-1.0'
#hf_seed = '5481190'
#hf_fmax = '10'
#hf_calpha = '-99'
#hf_accdir = 'acc'

# verify input incl. params.py
verify_binaries([hf_sim_bin])
verify_files([stat_file, hf_slip, hf_v_model])
verify_strings([hf_prefix, hf_t_len, hf_dt, hf_vs_moho, hf_fa_sig_1, hf_rv_sig_1, \
        hf_sdrop, hf_kappa, hf_qfexp, hf_rayset, hf_rvfac, hf_shal_rvfac, hf_deep_rvfac, \
        hf_czero, hf_site_amp, hf_mom, hf_rupv, hf_seed, hf_fmax, hf_calpha])
verify_user_dirs([hf_sim_dir])
verify_logfiles([local_statfile])
verify_user_dirs([hf_accdir],reset=True)

procs = 64
statfile_base = os.path.join(hf_sim_dir, 'local.tmp.sf_')
statfiles = []
for p in xrange(procs):
    statfiles.append('%s%d' % (statfile_base, p))

out_prefix = os.path.join(hf_accdir, hf_prefix)

# copy line not starting with '#' into a local stat_file copy
# file pointers for each processes' stat_file portion
fps = []
# number of stations for each process
nss = [0] * procs
for f in statfiles:
    fps.append(open(f, 'w'))
with open(stat_file, 'r') as fp:
    # in fp.readline() instead of in fp for AIX/old python? compatibility
    for line in fp.readlines():
        if line[0] != '#':
            fps[sum(nss) % procs].write(line)
            nss[sum(nss) % procs] += 1
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

set_permission(hf_sim_dir)

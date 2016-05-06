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

import os
from subprocess import call, Popen, PIPE
from shared import *
from params import *

local_statfile = os.path.join(hf_sim_dir, 'local.statfile')

# verify input incl. params.py
verify_binaries([hf_sim_bin])
verify_files([stat_file, hf_slip, hf_v_model])
verify_strings([hf_prefix, hf_t_len, hf_dt, hf_vs_moho, hf_fa_sig_1, hf_rv_sig_1, \
        hf_sdrop, hf_kappa, hf_qfexp, hf_rayset, hf_rvfac, hf_shal_rvfac, hf_deep_rvfac, \
        hf_czero, hf_site_amp, hf_mom, hf_rupv, hf_seed, hf_fmax, hf_calpha])
verify_dirs([hf_sim_dir])
verify_logfiles([local_statfile])
verify_user_dirs([hf_accdir])

out_prefix = os.path.join(hf_accdir, hf_prefix)

# copy line not starting with '#' into a local stat_file copy
n_stat = 0
with open(local_statfile, 'w') as lp:
    with open(stat_file, 'r') as fp:
        # in fp.readline() instead of in fp for AIX/old python? compatibility
        for line in fp.readlines():
            if line[0] != '#':
                lp.write(line)
                n_stat += 1
        # file must end with new line
        lp.write('\n')
# some info echoed in original script, might be useful for some
print("%d %s" % (n_stat, hf_v_model))

with open('hf.out', 'w') as op:
    # execute calculation, stdout and errout to hf.out file, stdin from pipe
    hf_pipe = Popen([hf_sim_bin], stdin = PIPE, stdout = op, stderr = op)
    # pipe input data into binary stdin
    # look at hf_sim_bin source code (same dir) for more info on parameters
    hf_pipe.communicate('\n%(hf_sdrop)s\n%(local_statfile)s\n%(out_prefix)s\n%(hf_rayset)s\n%(hf_site_amp)s\n4   0  0.02  19.9\n%(hf_seed)s\n%(n_stat)s\n%(hf_t_len)s %(hf_dt)s %(hf_fmax)s %(hf_kappa)s %(hf_qfexp)s\n%(hf_rvfac)s %(hf_shal_rvfac)s %(hf_deep_rvfac)s %(hf_czero)s %(hf_calpha)s\n%(hf_mom)s %(hf_rupv)s\n%(hf_slip)s\n%(hf_v_model)s\n%(hf_vs_moho)s\n-99 0.0 0.0 0.0 0.0 1\n-1\n%(hf_fa_sig_1)s 0.0 %(hf_rv_sig_1)s\n1\n' % locals())


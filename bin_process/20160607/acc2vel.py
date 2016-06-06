#!/usr/bin/env python2
"""
Converts HF acceleration seismograms to velocity versions.
TODO: is this file necessary? make generic? functionality included in match_seismo.

USAGE: Execute script from current directory.
    Directory should contain populated accelerations folder.

Converted to Python.
@date 06 May 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
"""
import sys
import os.path
sys.path.append(os.path.abspath(os.path.curdir))
from subprocess import call
from shared import *
from params import *

# verify input incl. params.py
verify_binaries([int_bin])
verify_files([stat_file])
verify_strings([hf_prefix])
verify_lists([int_comps])
verify_user_dirs([hf_accdir, hf_veldir])

stations = get_stations(stat_file)

for stat in stations:
    for comp in int_comps:
        file_in = os.path.join(hf_accdir, "%s_%s.%s" % (hf_prefix, stat, comp))
        file_out = os.path.join(hf_veldir, "%s.%s" % (stat, comp))
        # integrate ACC to get VEL
        call([int_bin, 'integ=1', "filein=%s" % file_in, \
                'inbin=0', 'outbin=0', "fileout=%s" % file_out])


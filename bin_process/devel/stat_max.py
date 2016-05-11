#!/usr/bin/env python2
"""
Generates LON LAT station list with maximum motion values.
Make sure the input station list contains the required stations.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 11 May 2016

USAGE: execute from current directory containing the 'Vel' folder.
    if outside $sim_dir, link within: 'ln -s location/to/stat_max.py simdir/'

ISSUES: 
"""

import os.path
from shared import *
from params import *

verify_files([FD_STATLIST])
verify_logfiles([MAX_STATFILE])
verify_strings([])
verify_user_dirs([])

stats, lats, lons = get_stations(FD_STATLIST, True)

# output
op = open(MAX_STATFILE, 'w')

for i, stat in enumerate(stats):
    velfile = os.path.join('Vel', '%s.000' % stat)
    if not os.path.isfile(velfile):
        print('WARNING: file not found: %s' % velfile)
        continue

    vp = open(velfile, 'r')
    # skip header
    vp.readline()
    vp.readline()

    max_val = 0.0
    for line in vp:
        for value in map(float, line.split()):
            if value > largest:
                largest = value

    op.write('%s %s %f' % (lons[i], lats[i], max_val))

op.close()


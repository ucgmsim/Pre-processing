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
from math import sqrt
from itertools import izip
from shared import *
from params import *

# return square sum
ss = lambda a, b: float(a) ** 2 + float(b) ** 2

verify_files([FD_STATLIST])
verify_logfiles([MAX_STATFILE])
verify_strings([])
verify_user_dirs([])

stats, lats, lons = get_stations(FD_STATLIST, True)
stat_count = len(stats)

# output
op = open(MAX_STATFILE, 'w')

for i, stat in enumerate(stats):
    # assume mostly no exceptions, use try block for speed
    try:
        v1p = open('Vel/%s.000' % stat, 'r')
        v2p = open('Vel/%s.090' % stat, 'r')
    except IOError:
        print('Could not open Vel files for station: %s' % stat)
        continue

    # skip headers, verification not necessary, speed more important
    v1p.readline()
    v1p.readline()
    v2p.readline()
    v2p.readline()

    max_value = 0.0
    cur_value = 0.0
    #TODO: which method is faster?
    # note assumption: files generated together so same no. of values per line
    for line in v1p:
        #for c000, c090 in izip(map(float, line.split()), map(float, v2p.readline().split())):
        cur_value = max(map(ss, line.split(), v2p.readline().split()))
        if cur_value > max_value:
            max_value = cur_value
            # no need to square root until finished
            #if c000**2 + c090**2 > max_val:
            #    max_val = c000**2 + c090**2

    v1p.close()
    v2p.close()
    op.write('%s %s %f\n' % (lons[i], lats[i], sqrt(max_value)))

op.close()
# new line after previous stdout.write()
print('')


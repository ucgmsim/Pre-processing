#!/usr/bin/env python2
"""
Generates LON,LAT,VEL IEEE754 32bit (single precision) binary file.
Make sure the input station list contains the required stations.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 June 2016

USAGE: 

ISSUES: 
"""

from math import sqrt

import h5py as h5
import numpy as np

# input
h5p = h5.File('virtual.hdf5', 'r')
num_stat = h5p.attrs['NSTAT']

# output
op = open('statgrid_max.bin', 'wb')

for stat_i, stat in enumerate(h5p):
    # show progress
    print('%d of %d' % (stat_i, num_stat))

    # creating variables is slow
    # '>f4' means big-endian, float, 32 bits (4 byte)
    np.array([ \
            h5p[stat].attrs['LON'], \
            h5p[stat].attrs['LAT'], \
            np.sqrt(np.max(np.add.reduce(np.power( \
                    h5p['%s/VEL' % (stat)][:,:2], 2), axis = 1))) \
            ], dtype='>f4').tofile(op)

op.close()
h5p.close()


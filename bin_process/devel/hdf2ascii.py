#!/usr/bin/env python2
"""
Generates ASCII VEL file for a gridpoint like winbin-aio but from HDF.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 June 2016

USAGE: first parameter is the X coordinate, second is Y.

ISSUES: 
"""

import sys

import h5py as h5
import numpy as np

# C code writes 6 values per line
VPL = 6

###
### INPUT verify command line arguments
###

try:
    x, y = map(int, sys.argv[1:3])
except ValueError:
    print 'Invalid input parameters, integers not found.'
    exit()

###
### INPUT verify file
###

h5p = h5.File('virtual.hdf5', 'r')
try:
    g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))
    group = h5p[g_name]
except KeyError:
    print('Gridpoint for X: %d, Y: %d not found. Exitting.' % (x, y))
    exit()

###
### INFORMATION which is in the ASCII header/rest of file
###

# start time? shift time? C code never modifies past initialising to 0
t_hr = 0
t_min = 0
# number of timesteps skipped times dt, -1 for starting at t = 0?
t_sec = -1.0
# epicentre distance, unset cmdline parameter by original winbin-aio
edist = 0.0
# following are never modified in C code past initialisation to 0
az = 0.0
baz = 0.0

# first entry in header
# leave as restricted short name (8 byte C-string) for compatibility
# put full name as title instead
name = group.attrs['NAME']
# entry not set by winbin-aio but structure allows it (normally blank)
title = g_name
# other values in second line
nt = h5p.attrs['NT']
dt = h5p.attrs['DT']

# load dataset as numpy array
data = h5p['%s/VEL' % (g_name)]
# dump to big endian binary file too
with open('%s.bin' % (name), 'wb') as bp:
    if sys.byteorder == 'little':
        data[...].byteswap().tofile(bp)
    else:
        data[...].tofile(bp)

###
### OUTPUT write to file
###

for comp_i in xrange(h5p.attrs['NCOMPS']):
    comp_name = h5p.attrs['COMP_%d' % (comp_i)]
    # header doesn't contain leading 0s
    try:
        comp_header = str(int(comp_name))
    except ValueError:
        comp_header = comp_name

    op = open('%s.%s' % (name, comp_name), 'w')
    # write header same format as C code
    op.write('%-10s %3s %s\n' % (name, comp_header, title))
    op.write('%d %12.5e %d %d %12.5e %12.5e %12.5e %12.5e\n' % \
            (nt, dt, t_hr, t_min, t_sec, edist, az, baz))
    # write values, also same format as C code
    # writing full lines turns a 7.4 second job into a 1.6 second job
    for full_line_i in xrange(nt // VPL):
        op.write('%13.5e%13.5e%13.5e%13.5e%13.5e%13.5e\n' % \
                tuple(data[full_line_i * VPL:full_line_i * VPL + VPL, comp_i]))
    # partial last line
    if nt % VPL != 0:
        for t in xrange(nt - nt % VPL, nt):
            op.write('%13.5e' % (data[t, comp_i]))
        op.write('\n')
    # close output
    op.close()

h5p.close()

#!/usr/bin/env python2
"""
Generates ASCII VEL file for a gridpoint like winbin-aio but from HDF.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 9 June 2016

USAGE: first parameter is the X coordinate, second is Y.

ISSUES: 
"""

import os
import sys
from zipfile import ZipFile

import h5py as h5
import numpy as np

from shared_bin import *
from tools import ll2gp, InputError

h5p = h5.File('virtual.hdf5', 'r')

###
### INPUT verify command line arguments
###

try:
    lat, lon = map(float, sys.argv[1:3])
except ValueError:
    print 'Invalid input parameters, float values not found.'
    exit()

###
### INPUT conversion lat lon to grid
###

# read required inputs
mlat = h5p.attrs['MLAT']
mlon = h5p.attrs['MLON']
nx = h5p.attrs['NX']
ny = h5p.attrs['NY']
dx = h5p.attrs['DX']
dy = h5p.attrs['DY']
hh = h5p.attrs['HH']
rot = h5p.attrs['ROT']

try:
    x, y = ll2gp(lat, lon, mlat, mlon, rot, nx, ny, hh, dx, dy)
except InputError:
    print('Input outside simulation domain.')
    exit()

###
### INPUT verify file
###

try:
    g_name = '%s%s' % (str(x).zfill(4), str(y).zfill(4))
    group = h5p[g_name]
except KeyError:
    print('Gridpoint for X: %d, Y: %d not found. Exitting.' % (x, y))
    exit()

###
### INFORMATION which is in the ASCII header/rest of file
###

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
    if NATIVE_ENDIAN == 'little':
        data[...].byteswap().tofile(bp)
    else:
        data[...].tofile(bp)

###
### OUTPUT write to file and ZIP
###

zf = ZipFile('%s.zip' % (name), 'w')
zf.write('%s/%s.bin' % ('.', name))
os.remove('%s/%s.bin' % ('.', name))
for comp_i in xrange(h5p.attrs['NCOMPS']):
    comp_name = h5p.attrs['COMP_%d' % (comp_i)]
    write_seis_ascii('.', data[:, comp_i], name, comp_name, nt, dt, title = title)
    zf.write('%s/%s.%s' % ('.', name, comp_name))
    os.remove('%s/%s.%s' % ('.', name, comp_name))
zf.close()
h5p.close()
# tell automation programs where the file is
print(zf.filename)

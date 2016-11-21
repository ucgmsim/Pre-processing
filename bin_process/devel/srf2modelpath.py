#!/usr/bin/env python2
"""
1.
Creates a binary file with LONG, LAT, VALUE.
Same format as a timeslice to be plotted in GMT.

2.
Creates a high resolution path of the fault
for use with GMT grdmask.

3.
Creates a standard corners file to plot fault planes.

ISSUES: crop will not work with multi-segment SRF
"""

from subprocess import call, Popen, PIPE

import numpy as np

from highres_modelpath import path_from_corners
from shared_srf import *
from tools import *

srf = 'h1.srf'

srf2xyz_bin = '/home/vap30/bin/srf2xyz'

# output data
out_file = 'slip_map.bin'
# values to extract (lon, lat, depth, value)
mask = np.array([True, True, False, True])

print("Loading SRF file %s (srf2xyz code)" % (srf))
proc = Popen([srf2xyz_bin, 'calc_xy=0', \
        'type=slip', 'nseg=-1', 'lonlatdep=1', \
        'dump_slip=0', 'infile=%s' % (srf)], \
        stdout = PIPE)
out, err = proc.communicate()
code = proc.wait()
# array will contain single values (should be 4)
array = np.fromstring(out, dtype = 'f4', sep = ' ')
# reshape to '4f' and remove depth
array = np.reshape(array, (len(array) // 4, 4))[:, mask]
print("Loading complete.")

###
### OUTPUT 1: binary file for GMT grid plotting
###
# store to file
array.astype(np.float32).tofile(out_file)

# find corners
corners = []
east, north = np.argmax(array[:, :2], axis = 0)
west, south = np.argmin(array[:, :2], axis = 0)
for extreme in [north, east, south, west]:
    corners.append(array[extreme, :2].tolist())

###
### OUTPUT 2: path around corners for GMT grdmask cropping
###
# list() creates a copy, prevent modification of local copy
path_from_corners(corners = list(corners), output = 'srf.modelpath_hr')

bounds = get_bounds(srf)

# find hypocentre
hypocentre = get_hypo(srf)

###
### OUTPUT 3: corners file for fault plane and hypocentre plot
###
# comment placed in corners file (above/below hypocentre)
# all lines (including the last one) must end with '\n'
corners_header = ('> header line here for specifics \n\
> This is the standard input file format where the \
hypocenter is first then for each \n\
>Hypocenter (reference??) \n', \
'> Below are the corners \
(first point repeated as fifth to close box \n')
with open('corners.txt', 'w') as cf:
    cf.write(corners_header[0])
    cf.write('%f %f\n' % (hypocentre[0], hypocentre[1]))
    cf.write(corners_header[1])
    for c, corners in enumerate(bounds):
        cf.write('> segment %d\n' % (c))
        for i in xrange(5):
            cf.write('%f %f\n' % tuple(corners[i % 4]))

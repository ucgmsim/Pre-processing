#!/usr/bin/env python2
"""
Creates a binary file with LONG, LAT, VALUE.
Same format as a timeslice to be plotted in GMT.

ISSUES: crop will not work with multi-segment SRF
"""

from subprocess import call, Popen, PIPE

import numpy as np

from highres_modelpath import path_from_corners

srf = '../../SrfGen/m7.80-240.0x77.0_s103246_v4.srf'

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

# store to file
array.astype(np.float32).tofile(out_file)

# find corners
corners = []
east, north = np.argmax(array[:, :2], axis = 0)
west, south = np.argmin(array[:, :2], axis = 0)
for extreme in [north, east, south, west]:
    corners.append(array[extreme, :2].tolist())

path_from_corners(corners = corners, output = 'srf.modelpath_hr')

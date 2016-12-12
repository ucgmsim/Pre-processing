#!/usr/bin/env python2
"""
1.
Creates binary files with LONG, LAT, VALUE for each segment.
Same format as a timeslice to be plotted in GMT.
Creates corresponding path files for masking area and GMT mask files.

2.
Creates a standard corners file to plot fault planes.

3.
Creates a GMT plot of the SRF file.
"""

import os
from shutil import rmtree

import numpy as np

from shared_gmt import *
from shared_srf import *
from tools import *

srf = 'v19.srf'

# output directory for srf resources
out_dir = os.path.abspath('SRF_RES')
if os.path.exists(out_dir):
    rmtree(out_dir)
os.makedirs(out_dir)

###
### OUTPUT 1: binary file for GMT grid plotting
###
print("Loading SRF file %s..." % (srf))
# get all corners
bounds = get_bounds(srf)
# TODO: fixup
np_bounds = np.array(bounds)
x_min, y_min = np.min(np.min(np_bounds, axis = 0), axis = 0)
x_max, y_max = np.max(np.max(np_bounds, axis = 0), axis = 0)
ll_region = (x_min, x_max, y_min, y_max)
for seg in xrange(len(bounds)):
    # create binary llv file for GMT overlay
    llslip = srf2llv(srf, seg = seg, type = 'slip', depth = False)
    llslip.astype(np.float32).tofile('%s/slip_map_%d.bin' % (out_dir, seg))
    # create a mask path for GMT overlay
    path_from_corners(corners = bounds[seg], min_edge_points = 100, \
            output = '%s/slip_map_%d.bounds' % (out_dir, seg))
    # create mask from path
    path = np.array(bounds[seg])
    x_min, y_min = np.min(path, axis = 0)
    x_max, y_max = np.max(path, axis = 0)
    grd_mask('%s/slip_map_%d.bounds' % (out_dir, seg), \
            '%s/slip_map_%d.mask' % (out_dir, seg), \
            dx = '0.1km', dy = '0.1km', \
            region = ll_region)
print("Loading complete.")


###
### OUTPUT 2: corners file for fault plane and hypocentre plot
###

# find hypocentre, use bounds from previous step
hypocentre = get_hypo(srf)

# comment placed in corners file (above/below hypocentre)
# all lines (including the last one) must end with '\n'
corners_header = ('> header line here for specifics \n\
> This is the standard input file format where the \
hypocenter is first then for each \n\
>Hypocenter (reference??) \n', \
'> Below are the corners \
(first point repeated as fifth to close box \n')
with open('%s/corners.txt' % (out_dir), 'w') as cf:
    cf.write(corners_header[0])
    cf.write('%f %f\n' % (hypocentre[0], hypocentre[1]))
    cf.write(corners_header[1])
    for c, corners in enumerate(bounds):
        cf.write('> segment %d\n' % (c))
        for i in xrange(5):
            cf.write('%f %f\n' % tuple(corners[i % 4]))


###
### OUTPUT 3: GMT MAP
###

makecpt('hot', '%s/slip.cpt' % (out_dir), 0, 150, 1, invert = True)
gmt_defaults(wd = out_dir)
p = GMTPlot('%s/plot.ps' % (out_dir))
p.background(11, 11)
p.spacial('M', ll_region, sizing = 6, left_margin = 1, bottom_margin = 2.5)
p.land()
p.water()
for seg in xrange(len(bounds)):
    p.overlay('%s/slip_map_%d.bin' % (out_dir, seg), \
            '%s/slip.cpt' % (out_dir), dx = '0.1km', dy = '0.1km', \
            climit = 2, crop_grd = '%s/slip_map_%s.mask' % (out_dir, seg), \
            land_crop = False)
p.finalise()
p.png()

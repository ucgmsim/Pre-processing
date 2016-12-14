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
import sys

import numpy as np

from shared_gmt import *
from shared_srf import *
from tools import *

srf = 'standard_m5.60-5.1x5.1_s103245.srf'
dx = '1k'
dy = '1k'
cpt = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
        'cpt', 'slip.cpt')
dpi = 300
if len(sys.argv) > 1:
    srf = sys.argv[1]
if len(sys.argv) > 3:
    dx, dy = sys.argv[2:4]

# output directory for srf resources
out_dir = os.path.abspath('SRF_RES')
if os.path.exists(out_dir):
    rmtree(out_dir)
os.makedirs(out_dir)

###
### OUTPUT 1: binary file for GMT grid plotting
###
print('Loading SRF file data from %s...' % (srf))
# get all corners
bounds = get_bounds(srf)
# gridding is computationally expensive
# each segment should include minimal region
seg_regions = []
# for percentile to automatically calculate colour palette range
values = np.array([], dtype = np.float32)
# find total extremes
np_bounds = np.array(bounds)
x_min, y_min = np.min(np.min(np_bounds, axis = 0), axis = 0)
x_max, y_max = np.max(np.max(np_bounds, axis = 0), axis = 0)
plot_region = (x_min - 0.1, x_max + 0.1, y_min - 0.1, y_max + 0.1)
# for plotting region on NZ-wide map
plot_bounds = '%f %f\n%f %f\n%f %f\n%f %f\n' % \
        (plot_region[0], plot_region[2], plot_region[1], plot_region[2], \
        plot_region[1], plot_region[3], plot_region[0], plot_region[3])
for seg in xrange(len(bounds)):
    # create binary llv file for GMT overlay
    llslip = srf2llv(srf, seg = seg, type = 'slip', depth = False)
    llslip.astype(np.float32).tofile('%s/slip_map_%d.bin' % (out_dir, seg))
    values = np.append(values, llslip[:, -1])
    # create a mask path for GMT overlay
    path_from_corners(corners = bounds[seg], min_edge_points = 100, \
            output = '%s/slip_map_%d.bounds' % (out_dir, seg))
    # create mask from path
    x_min, y_min = np.min(np_bounds[seg], axis = 0)
    x_max, y_max = np.max(np_bounds[seg], axis = 0)
    seg_regions.append((x_min, x_max, y_min, y_max))
    grd_mask('%s/slip_map_%d.bounds' % (out_dir, seg), \
            '%s/slip_map_%d.mask' % (out_dir, seg), \
            dx = dx, dy = dy, \
            region = seg_regions[seg])
percentile = np.percentile(values, 95)
print('Loading complete.')


###
### OUTPUT 2: corners file for fault plane and hypocentre plot
###
print('Creating corners file...')
# find hypocentre, use bounds from previous step
hypocentre = get_hypo(srf)
# standard corners file format
with open('%s/corners.txt' % (out_dir), 'w') as cf:
    cf.write('>hypocentre\n')
    cf.write('%f %f\n' % (hypocentre[0], hypocentre[1]))
    cf.write('>corners for each segment\n')
    for c, corners in enumerate(bounds):
        cf.write('>segment %d\n' % (c))
        for i in xrange(5):
            cf.write('%f %f\n' % tuple(corners[i % 4]))
print('Corners written to %s.' % ('%s/corners.txt' % (out_dir)))


###
### OUTPUT 3: GMT MAP
###
print('Plotting SRF on map...')
makecpt(cpt, '%s/slip.cpt' % (out_dir), 0, percentile, 1)
gmt_defaults(wd = out_dir)
p = GMTPlot('%s/srf_map.ps' % (out_dir))
p.background(20, 15)
p.spacial('M', plot_region, sizing = 6, left_margin = 1, bottom_margin = 2.5)
p.land()
p.water()
for seg in xrange(len(bounds)):
    p.overlay('%s/slip_map_%d.bin' % (out_dir, seg), \
            '%s/slip.cpt' % (out_dir), dx = dx, dy = dy, \
            climit = 2, crop_grd = '%s/slip_map_%s.mask' % (out_dir, seg), \
            land_crop = False, custom_region = seg_regions[seg])
p.fault(srf, is_srf = True)
p.coastlines()
p.sites(sites.keys())
p.text(plot_region[0], plot_region[3], os.path.basename(srf), \
        align = 'LB', size = '20p', dy = 0.1)
p.ticks(major = '10m', minor = '5m', sides = 'ws')

# draw NZ wide map to show rough location in country
p.spacial('M', nz_region, sizing = 4, left_margin = 7)
p.land()
p.water()
p.path(plot_bounds, is_file = False, close = True)
# try to make a line that extends from the box to the original map
# dots per degree longitude
dpdl = (4 * dpi) / (nz_region[1] - nz_region[0])
# gap between maps is 1 inch, dpi = pixel gap
# calculate equivalent longitude at first map
el = nz_region[0] - ((1.0 * dpi) / dpdl)
p.path('%f %f\n%f %f\n' % (plot_region[1], plot_region[2], el, nz_region[2]), \
        is_file = False)
p.fault(srf, is_srf = True, plane_width = '0.4p', top_width = '0.6p', \
        hyp_colour = 'red')
p.coastlines()
p.ticks(major = '2d', minor = '30m', sides = 'ws')

p.finalise()
p.png(dpi = dpi)
print('SRF plot at %s.' % ('%s/srf_map.png' % (out_dir)))

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

# can specify here or pass as command line argument
srf = 'standard_m5.60-5.1x5.1_s103245.srf'

dpi = 300

cpt = os.path.join(os.path.dirname(os.path.abspath(__file__)), \
        'cpt', 'slip.cpt')

if len(sys.argv) > 1:
    srf = sys.argv[1]

dx, dy = srf_dxy(srf)
dx = '%sk' % (dx)
dy = '%sk' % (dy)

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
maximum = np.max(values)
average = np.average(values)
subfaults = len(values)
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
# this is how hight the NZ map will end up being
# update if changed size from 4 inches wide or nz_region changed
p.spacial('M', nz_region, sizing = 4, \
        left_margin = 1, bottom_margin = 2.5)
# height of NZ map
full_height = mapproject(nz_region[0], nz_region[3], wd = out_dir)[1]
zoom_width = 6
# match height of zoomed in map with full size map
while True:
    p.spacial('M', plot_region, sizing = zoom_width)
    # height of this map
    zoom_height = mapproject(plot_region[0], plot_region[3], wd = out_dir)[1]
    if zoom_height > full_height * 1.01:
        zoom_width -= (zoom_height - full_height) / 2.0
    elif zoom_height < full_height * 0.99:
        zoom_width += (full_height - zoom_height) / 2.0
    else:
        break
p.land()
p.water()
for seg in xrange(len(bounds)):
    p.overlay('%s/slip_map_%d.bin' % (out_dir, seg), \
            '%s/slip.cpt' % (out_dir), dx = dx, dy = dy, \
            climit = 2, crop_grd = '%s/slip_map_%s.mask' % (out_dir, seg), \
            land_crop = False, custom_region = seg_regions[seg], \
            transparency = 30)
p.fault(srf, is_srf = True, hyp_colour = 'red')
p.coastlines()
p.sites(sites.keys())
# work out an ideal tick increment (ticks per inch)
# x axis is more constrainnig
major_tick = 0.05
while ((plot_region[1] - plot_region[0]) / major_tick) / zoom_width > 1.4:
    major_tick += 0.05
p.ticks(major = '%sd' % (major_tick), \
        minor = '%sd' % (major_tick / 2.0), sides = 'ws')

# draw NZ wide map to show rough location in country
p.spacial('M', nz_region, sizing = 4, left_margin = zoom_width + 1)
# height of NZ map
full_height = mapproject(nz_region[0], nz_region[3], wd = out_dir)[1]
p.land()
p.water()
p.path(plot_bounds, is_file = False, close = True, colour = 'blue')
# get displacement of box to draw zoom lines later
window_bottom = mapproject(plot_region[1], plot_region[2], wd = out_dir)
window_top = mapproject(plot_region[1], plot_region[3], wd = out_dir)
# also draw fault planes / hypocentre
p.fault(srf, is_srf = True, plane_width = '0.4p', top_width = '0.6p', \
        hyp_colour = 'red')
p.coastlines()
p.ticks(major = '2d', minor = '30m', sides = 'ws')

# draw zoom lines that extend from view box to original plot
p.spacial('X', \
        (0, window_bottom[0] + 1, 0, max(zoom_height, full_height)), \
        sizing = '%s/%s' % (window_top[0] + 1, max(zoom_height, full_height)), \
        left_margin = -1.0)
p.path('%f %f\n%f %f\n' % \
        (0, 0, window_bottom[0] + 1, window_bottom[1]), width = '0.6p', \
        is_file = False, split = '-', straight = True, colour = 'blue')
p.path('%f %f\n%f %f\n' % \
        (0, zoom_height, window_top[0] + 1, window_top[1]), width = '0.6p', \
        is_file = False, split = '-', straight = True, colour = 'blue')

# add text and colour palette
# position to enclose both plots
total_width = zoom_width + 1 + 4
total_height = max(zoom_height, full_height)
p.spacial('X', (0, total_width, 0, total_height + 2), \
        sizing = '%s/%s' % (total_width, total_height + 2), \
        left_margin = -zoom_width)
# SRF filename
p.text(total_width / 2.0, total_height, os.path.basename(srf), \
        align = 'CB', size = '20p', dy = 0.8)
# max slip
p.text(zoom_width / 2.0, total_height, 'Maximum slip: ', \
        align = 'RB', size = '14p', dy = 0.5)
p.text(zoom_width / 2.0 + 0.1, total_height, float('%.4f' % (maximum)), \
        align = 'LB', size = '14p', dy = 0.5)
# 95th percentile
p.text(zoom_width / 2.0, total_height, '95th percentile: ', \
        align = 'RB', size = '14p', dy = 0.3)
p.text(zoom_width / 2.0 + 0.1, total_height, float('%.4f' % (percentile)), \
        align = 'LB', size = '14p', dy = 0.3)
# average slip
p.text(zoom_width / 2.0, total_height, 'Average slip: ', \
        align = 'RB', size = '14p', dy = 0.1)
p.text(zoom_width / 2.0 + 0.1, total_height, float('%.4f' % (average)), \
        align = 'LB', size = '14p', dy = 0.1)
# planes
p.text(total_width - 4 / 2.0, total_height, 'Planes: ', \
        align = 'RB', size = '14p', dy = 0.5)
p.text(total_width - 4 / 2.0 + 0.1, total_height, len(bounds), \
        align = 'LB', size = '14p', dy = 0.5)
# dx and dy
p.text(total_width - 4 / 2.0, total_height, 'dX, dY: ', \
        align = 'RB', size = '14p', dy = 0.3)
p.text(total_width - 4 / 2.0 + 0.1, total_height, '%s, %s' % (dx, dy), \
        align = 'LB', size = '14p', dy = 0.3)
# subfaults
p.text(total_width - 4 / 2.0, total_height, 'Subfaults: ', \
        align = 'RB', size = '14p', dy = 0.1)
p.text(total_width - 4 / 2.0 + 0.1, total_height, subfaults, \
        align = 'LB', size = '14p', dy = 0.1)
# scale
p.cpt_scale(zoom_width / 2.0, -0.5, '%s/slip.cpt' % (out_dir), \
        int(percentile / 4.0), int(percentile / 4.0) / 2.0, \
        label = 'Slip (cm)', length = zoom_width)

p.finalise()
p.png(dpi = dpi)
print('SRF plot at %s.' % ('%s/srf_map.png' % (out_dir)))

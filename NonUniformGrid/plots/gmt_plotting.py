#!/usr/bin/env python2
"""
Script to plot the non-uniform grid with GMT
"""

import qcore_path
from qcore import gmt
import sys

# One way of getting the points on a file as they are needed
# awk '{print $1" "$2} input_file.ll > output_file.txt'

if len(sys.argv) < 2:
    print "Please provide a filename containing the Lon Lat points to plot"
    print "Usage: %s filename [corners_file] [output_name]" % sys.argv[0]
    exit(1)

name = "non_uniform_grid"
try:
    name = sys.argv[3]
except:
    print "Using default name %s" % name

p = gmt.GMTPlot(name+'.ps')
# background areas white, extra will be cropped
p.background(11, 11)
x_min = 166.2
x_max = 174.5
y_min = -47.4
y_max = -40.3

try:
    corners_file = sys.argv[2]
    with open(corners_file) as f:
        line = f.readlines()
        x_min, x_max, y_min, y_max = map(float,line[0].split())
        print "Using coordinates: %f %f %f %f" % (x_min, x_max, y_min, y_max)
except IndexError:
    print "Using default coordinates"

# mercator, 6 inches wide, leave space on page for tick labels
p.spacial('M', (x_min, x_max, y_min, y_max), sizing = 6, x_shift = 1, y_shift = 1)

# no topo
p.land(fill = 'darkgreen')
p.water()

# load the 1st argument containing all the points
p.points(sys.argv[1], is_file = True, shape = 'C', size = 0.01, fill = 'black', line = None)
# title
p.text((x_min + x_max) / 2, y_max, 'Non-Uniform Grid', dy = 0.2, align = 'CB', size = 20)

p.ticks(major = '1d', minor = '20m')

p.finalise()
p.png(dpi = 300, clip = True)

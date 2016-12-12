#!/usr/bin/env python2

import os
from shared_gmt import *

title = "title text"
x_min = 168.2
x_max = 177.9
y_min = -45.7
y_max = -37.85
dpi = 300
png_dir = os.path.abspath('Png')
topo_file = os.path.abspath('nztopo.grd')
gmt_temp = os.path.abspath('GMT_WD')

ll_region = (x_min, x_max, y_min, y_max)
ll_avg = sum(ll_region[:2]) / 2, sum(ll_region[2:]) / 2

cpt_land = '%s/land.cpt' % (gmt_temp)
template_bottom = '%s/bottom.ps' % (gmt_temp)

# topography colour scale
makecpt('%s/cpt/palm_springs_1.cpt' % (os.path.dirname(__file__)), cpt_land, \
        -250, 9000, inc = 10, invert = True)

b = GMTPlot(template_bottom)
# background can be larger as whitespace is later cropped
b.background(11, 11)
b.spacial('M', ll_region, sizing = 7, \
        left_margin = 1, bottom_margin = 2.5)
# title, fault model and velocity model subtitles
b.text(ll_avg[0], y_max, title, size = 20, dy = 0.2)
# topo, water, overlay cpt scale
b.topo(topo_file, cpt = cpt_land)
b.water(colour = 'lightblue', res = 'f')
# stations
b.points(os.path.abspath('xy.ll'), shape = 't', size = 0.08, \
        fill = 'white', line = 'white', line_thickness = 0.8)

t = b
t.coastlines()
t.seismo(os.path.abspath('gmt-seismo.xy'), 10000000, fmt = 'inc')
# simulation domain
#t.path(cnr_str, is_file = False, split = '-', \
#        close = True, width = '0.4p', colour = 'black')
# fault file - creating direct from SRF is slower
# OK if only done in template
#t.fault(srf_file, is_srf = True)
# ticks on top otherwise parts of map border may be drawn over
t.ticks(major = '60m', minor = '30m', sides = 'ws')

# render and do not continue
t.finalise()
t.png(dpi = dpi, out_dir = png_dir, clip = True)


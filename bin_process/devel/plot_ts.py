#!/usr/bin/env python2

# Created: 21 November 2016 from bash version (created 21 April 2016)
# Purpose: Generate visualisations of timeslices (PNG)
# Replacing Purpose 21-04-2016: Use bash instead of csh. Source vars from e3d.par
# Replacing Purpose 21-11-2016: Use python instead of bash for greater flexibility/integration.
# Authors: Viktor Polak <viktor.polak@canterbury.ac.nz>

# USAGE:
# Execute with python: "$ ./plot_ts.py" or "$ python2 plot_ts.py"
# Optional first parameter: specify number of threads
#   Default is user interactive with autodetect, capped to 8 threads

# ISSUES:
# could validate parameters/check if folders/files exist.

import multiprocessing as mp
import os
from shutil import rmtree
import sys

from params import *
from shared_gmt import *

# general setup
title = 'Mw7.8 14 Nov 2016 Kaikoura Earthquake'
fault_model = 'InSAR source model (modified from Hamling) v2'
vel_model = 'SIVM v1.65 h=0.4km'
# default is inches
map_width = 6
prefix = ts_out_prefix
ps_dir = os.path.join(sim_dir, 'PlotFiles')
png_dir = os.path.join(sim_dir, 'Png')
# temporary working directories for gmt are within here
# prevent multiprocessing issues by isolating processes
gmt_temp = os.path.join(sim_dir, 'GMT_WD')
dpi = 96
grd_dx = '1k'
grd_dy = '1k'

# default purpose is create multiple images
purpose = 'multi'
# threading: use first parameter if available (prevents user interaction)
if len(sys.argv) >= 2:
    try:
        processes = int(sys.argv[1])
        print('Using %s processes from cmdline.' % (processes))
    except ValueError:
        if sys.argv[1] == 'single':
            purpose = 'single'
            source = sys.argv[2]
            if not os.path.exists(source):
                print('Source: \'%s\' not found.' % (source))
                exit()
elif purpose == 'multi':
    virtual_cores = mp.cpu_count()
    if virtual_cores > 8:
        # performance issues with more threads (disk IO?)
        virtual_cores = 8

    print "Run on how many processes? [%d]: " % (virtual_cores),
    wanted_procs = raw_input()
    if wanted_procs == '':
        processes = virtual_cores
        print('Using default, %d processes.' % (virtual_cores))
    else:
        try:
            if int(wanted_procs) > 0:
                processes = int(wanted_procs)
                print('Running on %d processes.' % (processes))
        except ValueError:
            print('Invalid input. Exiting.')
            exit()


# definition of known/predetermined regions to plot for
# region to plot
# sites to display
# distance scale definition
if plot_region == 'CANTERBURY':
    x_min, x_max, y_min, y_max = 171.75, 173.00, -44.00, -43.20
    region_sites = ['Rolleston', 'Darfield', 'Lyttelton', 'Akaroa', \
            'Kaiapoi', 'Rakaia', 'Oxford']
    dist_scale = '-L172.50/-43.90/$modellat/25.0 -Ba30mf30mWSen'
elif plot_region == 'WIDERCANT':
    x_min, x_max, y_min, y_max = 170.52, 173.67, -44.4, -42.53
    region_sites = ['Rolleston', 'Darfield', 'Lyttelton', 'Akaroa', \
            'Kaiapoi', 'Rakaia', 'Oxford']
    dist_scale = '-L172.50/-43.90/$modellat/25.0 -Ba30mf30mWSen'
elif plot_region == 'SOUTHISLAND':
    x_min, x_max, y_min, y_max = 166.0, 174.5, -47.50, -40.00
    region_sites = ['Queenstown', 'Dunedin', 'Tekapo', 'Timaru', \
            'Christchurch', 'Haast', 'Greymouth', 'Westport', \
            'Kaikoura', 'Nelson', 'Blenheim']
    dist_scale = '-L173/-47/$modellat/50.0 -Ba60mf60mWSen'
elif plot_region == 'MIDNZ':
    x_min, x_max, y_min, y_max = 168.2, 177.9, -45.7, -37.85
    region_sites = ['Queenstown', 'Tekapo', 'Timaru', 'Christchurch', \
            'Haast', 'Greymouth', 'Westport', 'Kaikoura', 'Nelson', \
            'Blenheim', 'Wellington', 'Palmerston North', \
            'Masterton', 'Napier', 'New Plymouth', 'Taupo', 'Rotorua']
    dist_scale = '-L176/-45/$modellat/100.0 -Ba60mf60mWSen'
else:
    # TODO: optionally read custom region
    print('Plotting Region Not Understood')
    exit()

# avg lon/lat (midpoint of plotting region)
ll_avg = (x_min + x_max) / 2.0, (y_min + y_max) / 2.0
# combined region
ll_region = (x_min, x_max, y_min, y_max)

# clear output dirs
for out_dir in [ps_dir, png_dir, gmt_temp]:
    if os.path.isdir(out_dir):
        rmtree(out_dir)
    os.makedirs(out_dir)

# make temp file with simulation boundaries
# with -45 degree rotation:
#   c2
# c1  c3
#   c4
corners = []
with open(vel_mod_params, 'r') as vmpf:
    lines = vmpf.readlines()
    # make sure they are read in the correct order at efficiency cost
    for corner in ['c1=', 'c2=', 'c3=', 'c4=']:
        for line in lines:
            if corner in line:
                corners.append(map(float, line.split()[1:3]))



###
### PLOTTING STARTS HERE
###
######################################################

print('Creating resources...')
makecpt('cpt/palm_springs_1.cpt', '%s/land.cpt' % (gmt_temp), \
        -250, 9000, inc = 10, invert = True)
makecpt('hot', '%s/motion.cpt' % (gmt_temp), 0, 80, inc = 10, invert = True)
grd_mask('../sim.modelpath_hr', '%s/modelmask.grd' % (gmt_temp), \
        dx = grd_dx, dy = grd_dy, region = ll_region)

print('Creating bottom template...')
p = GMTPlot('%s/bottom.ps' % (gmt_temp))
p.background(11, 11)
p.spacial('M', ll_region, sizing = map_width, \
        left_margin = 1, bottom_margin = 2.5)
p.ticks(major = '60m', minor = '30m', sides = 'ws')
# title, fault model and velocity model subtitles
p.text(ll_avg[0], y_max, title, size = 20, dy = 0.6)
p.text(x_min, y_max, fault_model, size = 14, align = 'LB', dy = 0.3)
p.text(x_min, y_max, vel_model, size = 14, align = 'LB', dy = 0.1)
# water must go first due to GMT fails described in library
p.water(colour = 'lightblue', res = 'f')
p.topo('../nztopo.grd', cpt = 'land.cpt')
p.cpt_scale(3, -0.5, 'motion.cpt', 10, 10, label = 'ground motion (cm/s)')
# stations
p.points('../geonet_stations_20161119.ll', shape = 't', size = 0.08, \
        fill = None, line = 'white', line_thickness = 0.8)
p.leave()

print('Creating top template...')
p = GMTPlot('%s/top.ps' % (gmt_temp), append = True)
p.sites(region_sites)
p.coastlines()
p.leave()

copyfile('%s/bottom.ps' % (gmt_temp), '%s/combined.ps' % (gmt_temp))
with open('%s/combined.ps' % (gmt_temp), 'a') as c:
    with open('%s/top.ps' % (gmt_temp), 'r') as t:
        c.write(t.read())

p = GMTPlot('%s/combined.ps' % (gmt_temp), append = True)
p.finalise()
p.png(dpi = 96, clip = True)



#rmtree(gmt_temp)

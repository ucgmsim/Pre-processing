#!/usr/bin/env python2

"""
Created: 21 November 2016 from bash version (created 21 April 2016)
Purpose: Generate visualisations of timeslices (PNG)
Replacing Purpose 21-04-2016: Use bash instead of csh. Source vars from e3d.par
Replacing Purpose 21-11-2016: Use python instead of bash for greater flexibility.
Authors: Viktor Polak <viktor.polak@canterbury.ac.nz>

USAGE:
Execute with python: "$ ./plot_ts.py" or "$ python2 plot_ts.py"

USE CASES:
1. create timeslice png files
execute with 1st parameter being number of processes else interactive choice
2. create render with custom input same format as timeslice ('3f4')
execute with 1st parameter = single, second = source file

NOTES:
Processes seem to be IO limited (will see speedup with more processes with SSD).

ISSUES:
could validate parameters/check if folders/files exist.
"""

import multiprocessing as mp
import os
from shutil import copy, rmtree
import sys
from tempfile import mkdtemp
from time import time

from highres_modelpath import path_from_corners
try:
    from params import *
    prefix = os.path.abspath(ts_out_prefix)
    base_dir = os.path.abspath(sim_dir)
except ImportError:
    print('params.py not found. can only plot in single mode.')
    base_dir = os.path.abspath('.')
from shared_gmt import *

# general setup
title = 'Mw7.8 14 Nov 2016 Kaikoura Earthquake PGV'
fault_model = 'InSAR source model (modified from Hamling) v2'
vel_model = 'SIVM v1.65 h=0.4km'
# pick a predefined region or specify manual values
plot_region = 'MIDNZ'
x_min, x_max, y_min, y_max = 0, 0, 0, 0
dpi = 96
# overlay grid resolution 'k' for killometres, 'e' for metres.
grd_dx = '1k'
grd_dy = '1k'
# larger overlay grd limit improves speed, too large = bad results
convergence_limit = 2
# overlay cpt range min - max, with inc colour increment
cpt = 'hot'
cpt_min = 0
cpt_max = 80
cpt_inc = 10
cpt_legend = 'ground motion (cm/s)'
topo_file = os.path.abspath('nztopo.grd')
stat_file = os.path.abspath('geonet_stations_20161119.ll')
# values below lowcut or above highcut removed, if both then keep outside
# TODO: only lowcut implemented
cpt_lowcut = 2.0
cpt_highcut = None
# crop overlay to land
land_crop = False
map_width = 6
png_dir = os.path.join(base_dir, 'Png')
# seismogram overlay data, thickness, colour and format
seismo_ll = os.path.abspath('gmt-seismo.xy')
seismo_line = '0.8p'
seismo_colour = 'magenta'
seismo_fmt = 'time'
# temporary working directories for gmt are within here
# prevent multiprocessing issues by isolating processes
gmt_temp = os.path.join(base_dir, 'GMT_WD')

### END OF USER DEFINED INPUT

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
    dist_scale = '-L172.50/-43.90/$modellat/25.0'
elif plot_region == 'WIDERCANT':
    x_min, x_max, y_min, y_max = 170.52, 173.67, -44.4, -42.53
    region_sites = ['Rolleston', 'Darfield', 'Lyttelton', 'Akaroa', \
            'Kaiapoi', 'Rakaia', 'Oxford']
    dist_scale = '-L172.50/-43.90/$modellat/25.0'
elif plot_region == 'SOUTHISLAND':
    x_min, x_max, y_min, y_max = 166.0, 174.5, -47.50, -40.00
    region_sites = ['Queenstown', 'Dunedin', 'Tekapo', 'Timaru', \
            'Christchurch', 'Haast', 'Greymouth', 'Westport', \
            'Kaikoura', 'Nelson', 'Blenheim']
    dist_scale = '-L173/-47/$modellat/50.0'
elif plot_region == 'MIDNZ':
    x_min, x_max, y_min, y_max = 168.2, 177.9, -45.7, -37.85
    region_sites = ['Queenstown', 'Tekapo', 'Timaru', 'Christchurch', \
            'Haast', 'Greymouth', 'Westport', 'Kaikoura', 'Nelson', \
            'Blenheim', 'Wellington', 'Palmerston North', \
            'Masterton', 'Napier', 'New Plymouth', 'Taupo', 'Rotorua']
    dist_scale = '-L176/-45/$modellat/100.0'
else:
    # optionally allow custom region
    print('Plotting region not understood, using manual selection.')

# avg lon/lat (midpoint of plotting region)
ll_avg = (x_min + x_max) / 2.0, (y_min + y_max) / 2.0
# combined region
ll_region = (x_min, x_max, y_min, y_max)

# clear output dirs
for out_dir in [png_dir, gmt_temp]:
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
cnr_str = '\n'.join([' '.join(map(str, cnr)) for cnr in corners])



###
### PLOTTING STARTS HERE - TEMPLATE
###
######################################################

print('========== CREATING TEMPLATE ==========')
cpt_land = '%s/land.cpt' % (gmt_temp)
cpt_overlay = '%s/motion.cpt' % (gmt_temp)
mask = '%s/modelmask.grd'% (gmt_temp)
template_bottom = '%s/bottom.ps' % (gmt_temp)
template_top = '%s/top.ps' % (gmt_temp)

### create resources that are used throughout the process
t0 = time()
# topography colour scale
makecpt('%s/cpt/palm_springs_1.cpt' % (os.path.dirname(__file__)), cpt_land, \
        -250, 9000, inc = 10, invert = True)
# overlay colour scale
makecpt(cpt, cpt_overlay, cpt_min, cpt_max, \
        inc = cpt_inc, invert = True)
# simulation area mask
path_from_corners(corners = corners, min_edge_points = 100, \
        output = '%s/sim.modelpath_hr' % (gmt_temp))
grd_mask('%s/sim.modelpath_hr' % (gmt_temp), mask, \
        dx = grd_dx, dy = grd_dy, region = ll_region)
print('Created resources (%.2fs)' % (time() - t0))

### create a basemap template which all maps start with
t1 = time()
b = GMTPlot(template_bottom)
# background can be larger as whitespace is later cropped
b.background(11, 11)
b.spacial('M', ll_region, sizing = map_width, \
        left_margin = 1, bottom_margin = 2.5)
# title, fault model and velocity model subtitles
b.text(ll_avg[0], y_max, title, size = 20, dy = 0.6)
b.text(x_min, y_max, fault_model, size = 14, align = 'LB', dy = 0.3)
b.text(x_min, y_max, vel_model, size = 14, align = 'LB', dy = 0.1)
# topo, water, overlay cpt scale
b.topo(topo_file, cpt = cpt_land)
b.water(colour = 'lightblue', res = 'f')
b.cpt_scale(3, -0.5, cpt_overlay, cpt_inc, cpt_inc, label = cpt_legend)
# stations
b.points(stat_file, shape = 't', size = 0.08, \
        fill = None, line = 'white', line_thickness = 0.8)
# do not finalise, close file
b.leave()
print('Created bottom template (%.2fs)' % (time() - t1))

### add overlay if we are only plotting one input source
if purpose == 'single':
    b.enter()
    b.overlay(os.path.abspath(source), cpt_overlay, \
            dx = grd_dx, dy = grd_dy, climit = convergence_limit, \
            min_v = cpt_lowcut, crop_grd = mask, \
            contours = cpt_inc, land_crop = land_crop)
    b.leave()

### create map data which all maps will have on top
t2 = time()
if purpose == 'single':
    t = b
    t.enter()
else:
    t = GMTPlot(template_top, append = True)
t.sites(region_sites)
t.coastlines()
# simulation domain
t.path(cnr_str, is_file = False, split = '-', \
        close = True, width = '0.4p', colour = 'black')
# fault file - creating direct from SRF is slower
# OK if only done in template
t.fault(srf_file, is_srf = True)
# ticks on top otherwise parts of map border may be drawn over
t.ticks(major = '60m', minor = '30m', sides = 'ws')
t.leave()
print('Created top template (%.2fs)' % (time() - t2))

# render and do not continue
if purpose == 'single':
    t.enter()
    t.finalise()
    t.png(dpi = dpi, out_dir = png_dir, clip = True)
    exit()

### estimate time savings
t3 = time()
print('Time saved per timeslice: ~%.2fs' % (time() - t0))
print('Time saved over %d timeslices: ~%.2fs' \
        % (ts_total - ts_start, (ts_total - ts_start) * (time() - t0)))

print('========== TEMPLATE COMPLETE ==========')


###
### PLOTTING CONTINUES - TIME SLICE LOOPING
###
######################################################

def render_slice(n):
    t0 = time()

    # prepare resources in separate folder
    # prevents GMT IO errors on its conf/history files
    swd = mkdtemp(prefix = 'ts%.4dwd_' % (n), dir = gmt_temp)
    # name of slice postscript
    ps = '%s/ts%.4d.ps' % (swd, n)

    # copy GMT setup and basefile
    copy('%s/gmt.conf' % (gmt_temp), swd)
    copy('%s/gmt.history' % (gmt_temp), swd)
    copyfile(template_bottom, ps)
    s = GMTPlot(ps, append = True)
    # add timeslice timestamp
    s.text(x_max, y_max, 't=%.2fs' % (n * 0.1), \
            align = 'RB', size = '14p', dy = 0.1)

    # if we are plotting the total displacement
    abs_max('%s_ts%.4d.0' % (prefix, n), \
            '%s_ts%.4d.1' % (prefix, n), \
            '%s_ts%.4d.2' % (prefix, n), \
            '%s_ts%.4d.X' % (prefix, n), \
            native = False)
    s.overlay('%s_ts%.4d.X' % (prefix, n), cpt_overlay, \
            dx = grd_dx, dy = grd_dy, climit = convergence_limit, \
            min_v = cpt_lowcut, crop_grd = mask, \
            contours = cpt_inc, land_crop = land_crop)

    # append top file
    s.leave()
    with open(ps, 'a') as c:
        with open(template_top, 'r') as t:
            c.write(t.read())
    s.enter()

    # add seismograms if wanted
    if os.path.exists(seismo_ll):
        s.seismo(seismo_ll, n, fmt = seismo_fmt, \
                colour = seismo_colour, width = seismo_line)

    # create PNG
    s.finalise()
    s.png(dpi = dpi, out_dir = png_dir, clip = True)

    print('timeslice %d complete in %.2fs' % (n, time() - t0))

if purpose == 'multi':
    ts0 = time()
    #for i in xrange(ts_start, ts_total):
    #    render_slice(i)
    pool = mp.Pool(processes)
    pool.map(render_slice, xrange(ts_start, ts_total))
    print('FINISHED TIMESLICE SEGMENT IN %.2fs' % (time() - ts0))
    print('AVERAGE TIMESLICE TIME %.2fs' % \
            ((time() - ts0) / (ts_total - ts_start)))

# clear all working files
rmtree(gmt_temp)

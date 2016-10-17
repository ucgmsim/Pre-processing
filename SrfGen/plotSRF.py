#!/usr/bin/env python2

from math import ceil
from os import remove, path, makedirs
from subprocess import call, Popen, PIPE

import numpy as np

from setSrfParams import *

title = 'GP14.3'
out_dir = 'PlotFiles'
inputs = ['Srf/bev01_s103246_M6.66-6.91_FL9.85-12-20-14-7-11-8_FW18.45-19.67.srf']
# placed in outdir with '.ps' suffix
outputs = ['bev01']
cpt_dir = './srf_cpt'
srf2xyz_bin = '/home/vap30/bin/srf2xyz'

if not path.exists(out_dir):
    os.makedirs(out_dir)

plots = ['slip', 'trise', 'rake']
labels = ['Slip (cm)', 'Rise Time (s)', 'Rake (deg)']
x_label = 'along strike (km)'
y_label = 'W (km)'
# colour palettes for plot types
plot_cpt = ['%s/y2r2b.cpt' % (cpt_dir), \
        '%s/c2b2m.cpt' % (cpt_dir), \
        '%s/rain6.cpt' % (cpt_dir)]

# rounds plot maximums to nearest value
# also used for tick increments (1/5th rounded value)
near_v = [100.0, 0.5, 5.0]

# colour palette properties
min_v = [0.0, 0.0, -600.0]
# rounds plot maximums to nearest value
near_v = [100.0, 0.5, 5.0]
# set plot maximums manually (-1 = automatic scale)
max_v = [-1, -1, -1]
# colour change increment when using manual scale setting
# automatic scale will divide to produce 5 bins
# tick increments are 1/5th of rounded near_v value
inc_v = [20.0, 0.3, 10.0]
# CPT can create continuous gradient
continuous = [True, True, False]
# CPT colour increment
cpt_increment = [40.0, 0.6, 20.0]

# TODO: ?
psh = [-0.1, -0.1, -0.3]

# TODO: option flags
tinit = [True, False, False]
rmean = [False, False, True]
rakes = [False, False, True]
# precision (decimal places)
prec = [0, 1, 0]
slip_wgt = [False, True, False]
# display colour bar legend
color_bar = [True, True, False]
# whether rake is -180-180 (True), otherwise 0-360 (False)
flip_rake_angles = [False, False, True]

# TODO: NO IDEA
DX = [0.1, 0.1666667, 3.3333333, 2.0, 2.0]
DY = [0.1, 0.1666667, 2.5000000, 2.5, 2.5]

# TODO: take from setSrfParams
segments = []
seg_len = [13.0]
seg_wid = [9.0]

dx_rake = 0.67
dy_rake = 0.67
use_avg_rake = 0

# 4 for Mw~5, 8 for ~6, 12 for ~7, 40 for ~8
kminch = 4.0

ncol = 6
# tick increments
x_tick = 2.0
y_tick = 2.0
# TODO: shift plot?
yz = 7.5
xz = 1.5
x_shift = 3.2
# contour interval
contour_int = 1.0

xpmid = 4.25
# srf2xyz parameter (C boolean)
calc_xy = 0

grid_file = 'r.grd'


# set GMT properties
call(['gmtset', 'MAP_FRAME_WIDTH', '0.05', \
        'FONT_LABEL', '11', \
        'FORMAT_GEO_MAP', 'D', \
        'MAP_TICK_LENGTH_PRIMARY', '0.03', \
        'FORMAT_FLOAT_OUT', '%lg', \
        'PROJ_LENGTH_UNIT', 'inch', \
        'PS_PAGE_ORIENTATION', 'PORTRAIT'])

# store srf2xyz output at once instead of re-running
# does not use much RAM
xyz = []
for i, infile in enumerate(inputs):
    for xyz_type in ['slip', 'trise', 'rake', 'tinit']:
        proc = Popen([srf2xyz_bin, 'calc_xy=%d' % (calc_xy), \
                'type=%s' % (xyz_type), 'nseg=-1', \
                'dump_slip=1', 'infile=%s' % (infile)], \
                stdout = PIPE)
        out, err = proc.communicate()
        code = proc.wait()
        array = np.fromstring(out, dtype = 'f8', sep = ' ')
        xyz.append(np.reshape(array, (len(array) // 4, 4)))
xyz = np.array(xyz)

# store max and minimums for re-use
min_max = np.zeros(shape = (len(plots), 2), dtype = 'f')

for m, output in enumerate(outputs):

    ps_file = '%s/%s.ps' % (out_dir, output)
    if path.exists(ps_file):
        remove(ps_file)
    psf = open(ps_file, 'w')

    # start with blank page
    # use -F() to avoid input file piping
    call(['pstext', '-F()', '-R0/8.5/0/11', '-N', '-O', '-K', \
            '-X%f' % (xz), '-Y%f' % (yz)], stdout = psf)

    xsh = 0.0
    ysh = 0.0

    tailored_cpt = []
    ainc = []
    tick_inc = []

    for t, plot in enumerate(plots):
        tailored_cpt.append('temp%d.cpt' % (t))

        if flip_rake_angles[t]:
            min_max[t, 0] = round(np.min(np.where(xyz[t, :, 2] > 180.0, \
                    xyz[t, :, 2] - 360, xyz[t, :, 2])), prec[t])
            min_max[t, 1] = round(np.max(np.where(xyz[t, :, 2] > 180.0, \
                    xyz[t, :, 2] - 360, xyz[t, :, 2])), prec[t])
        else:
            min_max[t, 0] = round(np.min(xyz[t, :, 2]), prec[t])
            min_max[t, 1] = round(np.max(xyz[t, :, 2]), prec[t])

        # colour palette max is rounded up
        plot_max = ceil(min_max[t, 1] / near_v[t]) * near_v[t]
        # tick increment is always 5th of autosize
        tick_inc = plot_max * 0.2
        if max_v[t] == -1:
            plot_min = min_max[t, 0]
            plot_inc = plot_max * 0.2
        else:
            plot_max = max_v[t]
            plot_min = min_v[t]
            plot_inc = cpt_increment[t]
        ainc.append(plot_inc)

        # generate colour palettes to use
        cptp = open(tailored_cpt[t], 'w')
        mkcpt_p = Popen(['makecpt', '-C%s' % (plot_cpt[t]), \
                '-T%f/%f/%f' % (plot_min, plot_max, plot_inc), \
                '%s' % ('-Z' * continuous[t])], stdout = PIPE)
        out, err = mkcpt_p.communicate()
        code = mkcpt_p.wait()
        with open(tailored_cpt[t], 'w') as cptf:
            cptf.write(out)

    for i, infile in enumerate(inputs):
        x_inch = seg_len[i] / float(kminch)
        y_inch = seg_wid[i] / float(kminch)
        x_mid = x_inch + 0.3
        ys_len = y_inch * 0.8
        ys_mid = ys_len * 0.5

        scale = '-JX%f/-%f' % (x_inch, y_inch)
        region = '-R0/%f/0/%f' % (seg_len[i], seg_wid[i])

        title_p = Popen(['pstext', '-N', scale, region, \
                '-D0.0/0.2', '-K', '-O'], stdin = PIPE, stdout = psf)
        title_p.communicate( \
                '%f 0 20 0 1 2 %s' % (0.5 * seg_len[i], title))
        code = title_p.wait()

        for p, plot in enumerate(plots):
            if i < ncol:
                w = 'W'
            else:
                w = 'w'

            if not (i + 1) % ncol or \
                    (i == len(inputs) - 1 and p == len(plots) - 1):
                s = 'S'
            else:
                s = 's'

            val = xyz[p, :, 2]
            if flip_rake_angles[p]:
                val = np.where(val > 180, val - 360, val)
            if slip_wgt[p]:
                w = xyz[p, :, 3]
            else:
                w = np.ones(len(val))
            avg = round(np.sum(val * w) / np.sum(w), prec[p])
            mn, mx = min_max[p]


for cpt in tailored_cpt:
    remove(cpt)

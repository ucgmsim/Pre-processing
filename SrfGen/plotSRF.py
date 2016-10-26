#!/usr/bin/env python2

from math import ceil
from os import remove, path, makedirs
from subprocess import call, Popen, PIPE

import numpy as np

from setSrfParams import *

title = 'Multi Segment Testing'
out_dir = 'PlotFiles'
infile = 'Srf/bev01_s103245.srf'
# placed in outdir with '.ps' suffix
output = 'bev01'
cpt_dir = './srf_cpt'
srf2xyz_bin = '/home/vap30/bin/srf2xyz'
if not path.exists(out_dir):
    os.makedirs(out_dir)

# which plots to make (slip, trise, rake)
plots = ['slip', 'trise', 'rake']
labels = ['Slip (cm)', 'Rise Time (s)', 'Rake (deg)']
x_label = 'along strike (km)'
y_label = 'W (km)'

# rounds plot maximums to nearest value
# also used for tick increments (1/5th rounded value)
near_v = [50.0, 0.5, 5.0]
# colour palette properties
min_v = [0.0, 0.0, -3600]
# set plot maximums manually (-1 = automatic scale)
max_v = [-1, -1, 3600]
# colour change increment when using manual scale setting
# automatic scale will divide to produce 5 bins
# tick increments are 1/5th of rounded near_v value
# CPT colour increment
num_inc = [40.0, 0.6, 0.0]
tick_inc = [20, 0.3, 0.0]
# modified base colour palettes
plot_cpt = ['%s/y2r2b.cpt' % (cpt_dir), \
        '%s/c2b2m.cpt' % (cpt_dir), \
        '%s/rain6.cpt' % (cpt_dir)]
# whether CPT is continuous
continuous = [True, True, False]

# remove mean
rmean = [False, False, True]
# precision (decimal places)
prec = [0, 1, 0]
slip_wgt = [False, True, False]
# whether rake is -180-180 (True), otherwise 0-360 (False)
flip_rake_angles = [False, False, True]

dx_rake = 0.67
dy_rake = 0.67
use_avg_rake = False

# image tick increments
x_tick = 2.0
y_tick = 2.0

# contour interval
contour_int = 1.0

# properties of segments to put on page
if TYPE != 4:
    # not a multi segment SRF
    segments = [-1]
    seg_len = [FLEN]
    seg_wid = [FWID]
    # plotting grid resolution
    dy = [DWID]
    dx = [DLEN]
else:
    # multi segment
    segments = range(sum(M_NSEG))
    # TODO: retrieved by calling a single function
    seg_len = [M_FLEN[m][s] for m in xrange(len(M_FLEN)) \
            for s in xrange(len(M_FLEN[m]))]
    seg_wid = [M_FWID[m][s] for m in xrange(len(M_FWID)) \
            for s in xrange(len(M_FWID[m]))]
    # plotting grid resolution
    dy = [M_DWID[m][s] for m in xrange(len(M_DWID)) \
            for s in xrange(len(M_DWID[m]))]
    dx = [M_DLEN[m][s] for m in xrange(len(M_DLEN)) \
            for s in xrange(len(M_DLEN[m]))]

ncol = len(segments)
nrow = len(plots)
# set GMT properties
call(['gmtset', 'MAP_FRAME_WIDTH', '0.025', \
        'FONT_LABEL', '11', \
        'FORMAT_GEO_MAP', 'D', \
        'MAP_TICK_LENGTH_PRIMARY', '0.06', \
        'FORMAT_FLOAT_OUT', '%lg', \
        'PROJ_LENGTH_UNIT', 'inch', \
        'PS_MEDIA', '=', 'A4', \
        'PS_PAGE_ORIENTATION', 'LANDSCAPE'])
# width x height of PS_MEDIA defined above
media = [11.7, 8.3]
# top space for title
top_margin = 1.0
# bottom space for tick labels
bottom_margin = 0.75
# left space for scale, margin
left_margin = 0.8
# right margin for colour scales
right_margin = 1.0
# vertical space between plots (leave room for headings)
row_spacing = 0.4
# horizontal space between plots
col_spacing = 0.2

# max Y axis plot space per plot after all margins
plot_y = (media[1] - (nrow - 1) * row_spacing \
        - top_margin - bottom_margin) / nrow
# kilometres per unit (inch) of paper
# based on Y axis limitations
km_inch = max(seg_wid) / plot_y
# required X axis space considering above
plot_x_total = left_margin + (ncol - 1) * col_spacing \
        + sum(seg_len) / km_inch + right_margin
# check if enough X axis space
if plot_x_total > media[0]:
    # max X axis plot space
    plot_x = media[0] - plot_x_total + sum(seg_len) / km_inch
    # recalculate based on more limited X axis limitations
    km_inch = sum(seg_len) / plot_x
# space to bottom of first plot (placed at top)
yz = media[1] - top_margin - plot_y

# store srf2xyz output at once instead of re-running
# does not use much RAM
print("Loading SRF file %s (srf2xyz code)" % (infile))
xyz = []
for s in segments:
    seg = []
    for xyz_type in ['slip', 'trise', 'rake', 'tinit']:
        proc = Popen([srf2xyz_bin, 'calc_xy=0', \
                'type=%s' % (xyz_type), 'nseg=%i' % (s), \
                'dump_slip=1', 'infile=%s' % (infile)], \
                stdout = PIPE)
        out, err = proc.communicate()
        code = proc.wait()
        array = np.fromstring(out, dtype = 'f4', sep = ' ')
        seg.append(np.reshape(array, (len(array) // 4, 4)))
    xyz.append(np.array(seg))
print("Loading complete.")

ps_file = '%s/%s.ps' % (out_dir, output)
if path.exists(ps_file):
    remove(ps_file)
psf = open(ps_file, 'w')

# initialise postscript
# place origin for first plot
title_p = Popen(['pstext', '-N', '-K', \
        '-D0.0/%s' % (plot_y + row_spacing), \
        '-JX%s/%s' % (media[0], media[1]), \
        '-R0/%s/0/%s' % (media[0], media[1]), \
        '-X%f' % (left_margin), '-Y%f' % (yz)], \
        stdin = PIPE, stdout = psf)
title_p.communicate( \
        '%f 0 20 0 1 2 %s' % (0.5 * seg_len[0], title))
title_p.wait()

xsh = 0.0
ysh = 0.0

tailored_cpt = []
# store overall max and minimums for re-use
min_max = np.zeros(shape = (len(plots), 2), dtype = 'f')
for p, plot in enumerate(plots):
    tailored_cpt.append('temp%d.cpt' % (p))

    # find range of values, store for later use
    if flip_rake_angles[p]:
        seg_min = [np.min(np.where(xyz[s][p, :, 2] > 180.0, \
                        xyz[s][p, :, 2] - 360, xyz[s][p, :, 2])) \
                                for s in xrange(len(segments))]
        seg_min = [np.max(np.where(xyz[s][p, :, 2] > 180.0, \
                        xyz[s][p, :, 2] - 360, xyz[s][p, :, 2])) \
                                for s in xrange(len(segments))]
    else:
        seg_min = [np.min(xyz[s][p, :, 2]) for s in xrange(len(segments))]
        seg_max = [np.max(xyz[s][p, :, 2]) for s in xrange(len(segments))]
    min_max[p, 0] = round(min(seg_min), prec[p])
    min_max[p, 1] = round(max(seg_max), prec[p])

    # automatically scale colour palette
    if max_v[p] == -1:
        max_v[p] = ceil(min_max[p, 1] / near_v[p]) * near_v[p]
        min_v[p] = min_max[p, 0]
        tick_inc[p] = max_v[p] * 0.1
        num_inc[p] = max_v[p] * 0.2

    # generate colour palettes to use
    cmd = ['makecpt', '-C%s' % (plot_cpt[p]), \
            '-T%s/%s/%s' % (min_v[p], max_v[p], tick_inc)]
    if continuous[p]:
        cmd.extend(['-Z'])
    mkcpt_p = Popen(cmd, stdout = PIPE)
    out, err = mkcpt_p.communicate()
    code = mkcpt_p.wait()
    with open(tailored_cpt[p], 'w') as cptf:
        cptf.write(out)

for s, seg in enumerate(segments):
    x_inch = seg_len[seg] / float(km_inch)
    y_inch = seg_wid[seg] / float(km_inch)
    x_mid = x_inch + 0.3
    ys_len = y_inch * 0.8
    ys_mid = ys_len * 0.5

    scale = '-JX%f/-%f' % (x_inch, y_inch)
    region = '-R0/%f/0/%f' % (seg_len[seg], seg_wid[seg])

    for p, plot in enumerate(plots):
        # which sides to put tick labels on (captial)
        if s == 0:
            W = 'W'
        else:
            W = 'w'
        if p == nrow - 1:
            S = 'S'
        else:
            S = 's'

        val = xyz[seg][p, :, 2]
        if flip_rake_angles[p]:
            val = np.where(val > 180, val - 360, val)
        if slip_wgt[p]:
            w = xyz[seg][p, :, 3]
        else:
            w = np.ones(len(val))
        avg = round(np.sum(val * w) / np.sum(w), prec[p])
        mn, mx = min_max[p]

        stuff = np.copy(xyz[seg][p, :, :3])
        if rmean[p]:
            stuff[:, 2] = val - avg
        else:
            stuff[:, 2] = val
        grdin = ''.join(['%13.5e %13.5e %13.5e\n' \
                % (stuff[t][0], stuff[t][1], stuff[t][2]) \
                for t in xrange(len(stuff))])
        grd_p = Popen(['xyz2grd', '-Ggrd.grd', '-I%f/%f' \
                % (dx[s], dy[s]), region, '-r'], stdin = PIPE)
        grd_p.communicate(grdin)
        code = grd_p.wait()

        # plot image / background
        call(['grdimage', 'grd.grd', scale, region, '-C%s' % tailored_cpt[p], \
                '-B%f:%s:/%f:%s:%s%sen' % (x_tick, x_label, y_tick, y_label, W, S), \
                '-K', '-O', '-X%f' % (xsh), '-Y%s' % (ysh)], stdout = psf)


        # add label
        lbl_p = Popen(['pstext', scale, region, '-N', '-O', '-K', \
                '-D0.025/0.05'], stdin = PIPE, stdout = psf)
        lbl_p.communicate('0.0 0.0 12 0 1 1 %s\n%s 0.0 11 0 0 3 %s / %s / %s\n' \
                % (labels[p], seg_len[seg], min_max[p][0], avg, min_max[p][1]))
        code = lbl_p.wait()


        # add scale
        if plots[p] != 'rake' and s == ncol - 1:
            call(['psscale', '-Ef', '-C%s' % (tailored_cpt[p]), \
                    '-D%s/%s/%s/0.1' % (x_mid, ys_mid, ys_len), \
                    '-K', '-O', '-Ba%sf%sg%s::/::' \
                    % (num_inc[p], tick_inc[p], tick_inc[p])], stdout = psf)


        # show TINIT
        if plots[p] == 'slip':
            data = ''.join(['%s %s %s\n' \
                    % (xyz[seg][3, t, 0], xyz[seg][3, t, 1], xyz[seg][3, t, 2]) \
                    for t in xrange(len(xyz[seg][3]))])
            grd_p = Popen(['xyz2grd', '-Ggrd.grd', '-I%f/%f' \
                    % (dx[s], dy[s]), region, '-r'], stdin = PIPE)
            grd_p.communicate(data)
            code = grd_p.wait()
            call(['grdcontour', '-A%s+gwhite@50' % contour_int, \
                    'grd.grd', scale, region, '-W1.0', '-K', '-O'], \
                    stdout = psf)

        # show RAKE
        if plots[p] == 'rake':
            data = xyz[seg][2]
            ## only show arrows at dx, dy resolution
            nx = int(ceil(seg_len[seg] / dx_rake))
            ny = int(ceil(seg_wid[seg] / dy_rake))
            gridpoints = nx * ny

            # x and y position on plot
            x0 = np.zeros(gridpoints)
            y0 = np.zeros(gridpoints)
            # rake angle and slip amount
            rk = np.zeros(gridpoints)
            sp = np.zeros(gridpoints)

            # decimated x, y and grid array positons
            ix = (data[:, 0] / dx_rake).astype('i4')
            iy = (data[:, 1] / dy_rake).astype('i4')
            ip = ix + iy * nx
            # use closest position to middle of grid square
            if not use_avg_rake:
                mr = ((ix + 0.5) * dx_rake - data[:, 0]) ** 2 \
                        + ((iy + 0.5) * dy_rake - data[:, 1]) ** 2
                for x in xrange(gridpoints):
                    mid = np.where(ip == x)[0][np.argmin(mr[ip == x])]
                    x0[x] = (ix[mid] + 0.5) * dx_rake
                    y0[x] = (iy[mid] + 0.5) * dy_rake
                    rk[x] = data[mid, 2]
                    sp[x] = data[mid, 3]
            # use average value of all positions in grid square
            else:
                x0[ip] = (ix + 0.5) * dx_rake
                y0[ip] = (iy + 0.5) * dy_rake
                for x in xrange(gridpoints):
                    rk[x] = np.average(data[ip == x, 2])
                    sp[x] = np.average(data[ip == x, 3])

            output = []
            for x in xrange(gridpoints):
                if not np.isnan(rk[x]):
                    output.append('%13.5e %13.5e %13.5e %f\n' % (x0[x], y0[x], \
                            -rk[x], 0.4 * sp[x] / min_max[0, 1]))
            output = ''.join(output)
            rake_p = Popen(['psxy', '-Sv0.005i/0.04i/0.02i', '-W2', '-G0/0/0', \
                    '-K', '-O', scale, region], stdin = PIPE, stdout = psf)
            out, err = rake_p.communicate(output)
            code = rake_p.wait()

        # draw plots on top of eachother
        xsh = 0
        ysh = -(y_inch + row_spacing)
    # move over for next segment
    xsh = x_inch + col_spacing
    # move to top for next segment
    ysh = (y_inch + row_spacing) * (nrow - 1)
psf.close()

for cpt in tailored_cpt:
    remove(cpt)

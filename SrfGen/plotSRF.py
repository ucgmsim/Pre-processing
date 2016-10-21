#!/usr/bin/env python2

from math import ceil
from os import remove, path, makedirs
from subprocess import call, Popen, PIPE

import numpy as np

from setSrfParams import *

title = 'GP14.3'
out_dir = 'PlotFiles'
inputs = ['m6.20-16.0x9.0_s1129571.srf']
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
max_v = [200, 3.0, 600]
# colour change increment when using manual scale setting
# automatic scale will divide to produce 5 bins
# tick increments are 1/5th of rounded near_v value
# CPT colour increment
num_inc = [40.0, 0.6, 20.0]
tick_inc = [20, 0.3, 10.0]
# whether CPT should be continuous (not discretised)
continuous = [True, True, False]

# TODO: ?
psh = [-0.1, -0.1, -0.3]

# remove mean
rmean = [False, False, True]
# precision (decimal places)
prec = [0, 1, 0]
slip_wgt = [False, True, False]
# display colour bar legend
colour_bar = [True, True, False]
# whether rake is -180-180 (True), otherwise 0-360 (False)
flip_rake_angles = [False, False, True]

# TODO: plotting grid resolution?
DX = [0.1, 0.1666667, 3.3333333, 2.0, 2.0]
DY = [0.1, 0.1666667, 2.5000000, 2.5, 2.5]

# TODO: take from setSrfParams
segments = []
seg_len = [16.0]
seg_wid = [9.0]

dx_rake = 0.67
dy_rake = 0.67
use_avg_rake = False

# 4 for Mw~5, 8 for ~6, 12 for ~7, 40 for ~8
# kilometres per inch of paper
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
        array = np.fromstring(out, dtype = 'f4', sep = ' ')
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
    call(['psxy', '-T', '-JX8.5/11.0', '-R0/8.5/0/11', \
            '-N', '-K', \
            '-X%f' % (xz), '-Y%f' % (yz)], stdout = psf)

    xsh = 0.0
    ysh = 0.0

    tailored_cpt = []
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

        # automatically scale colour palette
        if max_v[t] == -1:
            max_v[t] = ceil(min_max[t, 1] / near_v[t]) * near_v[t]
            min_v[t] = min_max[t, 0]
            tick_inc[t] = max_v[t] * 0.1
            num_inc[t] = max_v[t] * 0.2

        # generate colour palettes to use
        cmd = ['makecpt', '-C%s' % (plot_cpt[t]), \
                '-T%s/%s/%s' % (min_v[t], max_v[t], tick_inc)]
        if continuous[t]:
            cmd.extend(['-Z'])
        mkcpt_p = Popen(cmd, stdout = PIPE)
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
                W = 'W'
            else:
                W = 'w'

            if not (i + 1) % ncol or \
                    (i == len(inputs) - 1 and p == len(plots) - 1):
                S = 'S'
            else:
                S = 's'

            val = xyz[p, :, 2]
            if flip_rake_angles[p]:
                val = np.where(val > 180, val - 360, val)
            if slip_wgt[p]:
                w = xyz[p, :, 3]
            else:
                w = np.ones(len(val))
            avg = round(np.sum(val * w) / np.sum(w), prec[p])
            mn, mx = min_max[p]

            val = xyz[p, :, 2]
            if flip_rake_angles[p]:
                val = np.where(val > 180, val - 360, val)
            stuff = np.copy(xyz[p, :, :3])
            if rmean[p]:
                stuff[:, 2] = val - avg
            else:
                stuff[:, 2] = val
            grdin = ''.join(['%13.5e %13.5e %13.5e\n' \
                    % (stuff[t][0], stuff[t][1], stuff[t][2]) \
                    for t in xrange(len(stuff))])
            grd_p = Popen(['xyz2grd', '-Ggrd.grd', '-I%f/%f' \
                    % (DX[i], DY[i]), region, '-r'], stdin = PIPE)
            grd_p.communicate(grdin)
            code = grd_p.wait()

            # plot image / background
            call(['grdimage', 'grd.grd', scale, region, '-C%s' % tailored_cpt[p], \
                    '-B%f:"%s":/%f:"%s":%s%sen' % (x_tick, x_label, y_tick, y_label, W, S), \
                    '-K', '-O', '-X%f' % (xsh), '-Y%s' % (ysh)], stdout = psf)


            # add label
            lbl_p = Popen(['pstext', scale, region, '-N', '-O', '-K', \
                    '-D0.025/0.05'], stdin = PIPE, stdout = psf)
            lbl_p.communicate('0.0 0.0 12 0 1 1 %s\n%s 0.0 11 0 0 3 %s / %s / %s\n' \
                    % (labels[p], seg_len[i], min_max[p][0], avg, min_max[p][1]))
            code = lbl_p.wait()


            # add scale
            if colour_bar[p]:
                call(['psscale', '-Ef', '-C%s' % (tailored_cpt[p]), \
                        '-D%s/%s/%s/0.1' % (x_mid, ys_mid, ys_len), \
                        '-K', '-O', '-Ba%sf%sg%s::/::' \
                        % (num_inc[p], tick_inc[p], tick_inc[p])], stdout = psf)


            # show TINIT
            if plots[p] == 'slip':
                data = ''.join(['%s %s %s\n' \
                        % (xyz[3, t, 0], xyz[3, t, 1], xyz[3, t, 2]) \
                        for t in xrange(len(xyz[3]))])
                grd_p = Popen(['xyz2grd', '-Ggrd.grd', '-I%f/%f' \
                        % (DX[i], DY[i]), region, '-r'], stdin = PIPE)
                grd_p.communicate(data)
                code = grd_p.wait()
                call(['grdcontour', 'grd.grd', scale, region, \
                        '-C%s' % contour_int, '-W1.0', '-K', '-O'], \
                        stdout = psf)

            # show RAKE
            if plots[p] == 'rake':
                data = xyz[2]
                ## only show arrows at dx, dy resolution
                nx = int(ceil(seg_len[i] / dx_rake))
                ny = int(ceil(seg_wid[i] / dy_rake))
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
                        '-K', '-O', '-J', '-R'], stdin = PIPE, stdout = psf)
                out, err = rake_p.communicate(output)
                code = rake_p.wait()

            xsh = 0.0
            if i + 1 % ncol == 0:
                xsh = x_inch + 0.2

            ysh = -(y_inch + 0.4)
            if i + 1 % ncol == 0:
                ysh = -ysh * (ncol - 1)
psf.close()

for cpt in tailored_cpt:
    remove(cpt)

#!/usr/bin/env python

from os import path, remove
from decimal import Decimal
from subprocess import call
from shared import *
from params import *

combined_stat_temp = 'all_stats.ll.tmp'
statgrid_temp = 'statgrid.tmp'
pre_ll_temp = 'll_pre.tmp'
post_ll_temp = 'll_post.tmp'
tmp_files = [combined_stat_temp, statgrid_temp, pre_ll_temp, post_ll_temp]

verify_binaries([LL2GRID_BIN, XYZ2LL_BIN])
verify_strings([MODEL_LAT, MODEL_LON, MODEL_ROT, hh, nx, ny, nz, X_BND, Y_BND, \
        CENTER_ORIGIN, GEOPROJ])
verify_lists([STAT_FILES])
verify_logfiles([combined_stat_temp, statgrid_temp, pre_ll_temp])
verify_dirs([stat_dir])

nx, ny = map(int, [nx, ny])

# default - dont change - why?
# seems to allow decimal gridpoints
VARGRID = 0

# define the length of the model domain in the local X and Y axes
XLEN = float(nx) * float(hh)
YLEN = float(ny) * float(hh)

# define the location of the Binary files, the output STATCORDS file (for the FD sim), an 'echo' file to
# check the input is read correctly (FDLLCORDS), and the name of the input file (STATFILES)
STATCORDS = path.join(stat_dir, 'fd_nz01-h%s.statcords' % (hh))
FDLLCORDS = path.join(stat_dir, 'fd_nz01-h%s.ll' % (hh))

with open(combined_stat_temp, 'w') as csp:
    for sf in STAT_FILES:
        with open(sf, 'r') as sp:
            csp.write(''.join(sp.readlines()))

with open ('duplicate_stats', 'w') as dp:
    call([LL2GRID_BIN, 'mlat=%s' % MODEL_LAT, 'mlon=%s' % MODEL_LON, 'rotate=%s' % MODEL_ROT, \
            'center_origin=%s' % CENTER_ORIGIN, 'geoproj=%s' % GEOPROJ, \
            'var_grid=%s' % VARGRID, 'xlen=%s' % XLEN, 'ylen=%s' % YLEN, \
            'nx=%s' % nx, 'ny=%s' % ny, 'h=%s' % hh, 'xbnd=%s' % X_BND, 'ybnd=%s' % Y_BND, \
            'infile=%s' % combined_stat_temp, 'outfile=%s' % statgrid_temp], \
            stdout = dp, stderr = dp)

if VARGRID == 1:
    # XMIN, XMAX, YMIN, YMAX variables were missing
    # X_BND, XLEN-X_BND, Y_BND, YLEN-XBND ??
    with open(statgrid_temp, 'r') as sp:
        # skip item count
        sp.readline()
        # output which should be printed after number of values
        out_lines = []
        for line in sp.readlines():
            try:
                x, y, z = map(float, line.split()[:3])
                stat = line.split()[3]
            except ValueError:
                continue
            if XMIN < x < XMAX and YMIN < y < YMAX:
                out_lines.append('%10.5f %10.5f %10.5f %s\n' % (x, y, z, stat))
    with open(STATCORDS, 'w') as scp:
        scp.write('%d\n' % len(out_lines))
        scp.write(''.join(out_lines))

else:
    with open(statgrid_temp, 'r') as sp:
        # skip item count
        sp.readline()
        # output which should be printed after number of values
        out_lines_statcords = []
        out_lines_temp = []
        for line in sp.readlines():
            try:
                x, y, z = map(int, line.split()[:3])
                stat = line.split()[3]
            except ValueError:
                continue
            if X_BND < x < nx - X_BND and Y_BND < y < ny - Y_BND:
                out_lines_statcords.append('%5d %5d %5d %s\n' % (x, y, z, stat))
                out_lines_temp.append('%10.4f %10.4f %10.4f %s\n' % \
                        (x * Decimal(hh), y * Decimal(hh), (z - 1) * Decimal(hh), stat))
    with open(STATCORDS, 'w') as scp:
        scp.write('%d\n' % len(out_lines_statcords))
        scp.write(''.join(out_lines_statcords))
    with open(pre_ll_temp, 'w') as pp:
        pp.write('%d\n' % len(out_lines_temp))
        pp.write(''.join(out_lines_temp))

call([XYZ2LL_BIN, 'mlat=%s' % MODEL_LAT, 'mlon=%s' % MODEL_LON, 'mrot=%s' % MODEL_ROT, \
        'center_origin=%s' % CENTER_ORIGIN, 'geoproj=%s' % GEOPROJ, \
        'xlen=%s' % XLEN, 'ylen=%s' % YLEN, 'infile=%s' % pre_ll_temp, 'outfile=%s' % post_ll_temp])

with open(post_ll_temp, 'r') as pp:
    with open(FDLLCORDS, 'w') as lp:
        for line in pp.readlines():
            items = line.split()
            lp.write('%11.5f %11.5f %s\n' % (float(items[0]), float(items[1]), items[3]))

for tmp in tmp_files:
    remove(tmp)

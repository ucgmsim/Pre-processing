#!/usr/bin/env python2
"""
Takes a model coords file and makes a statgrid file.
Statgrid files can be used for creating shakemap type plots of
the max shaking over the full ground motion time series.

USAGE: Execute in directory with shared.py and params.py linked in.
    Remember to edit params.py as required.

Converted to Python. CSH ~ 5.15 seconds vs Python ~ 5.40 seconds.
@date 10 May 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
"""
import sys
import os.path
from shared import *

sys.path.append(os.path.abspath(os.path.curdir))
from params import *

# for testing
#MODEL_COORDS = 'model_coords_nz01-h0.100'
verify_files([MODEL_COORDS])
verify_logfiles(['%s.ll' % (STATGRID_GEN), '%s.statcords' % (STATGRID_GEN)])
verify_strings([nx, ny, dx_ts, dy_ts, X_BND_PAD, Y_BND_PAD])

nx, ny, dx_ts, dy_ts = map(int,[nx, ny, dx_ts, dy_ts])

llp = open('%s.ll' % (STATGRID_GEN), 'w')
sc = []
with open(MODEL_COORDS, 'r') as mcp:
    coords = mcp.readlines()

for line in coords:
    info = line.split()
    # it is faster to calculate each time than creating separate variables
    # lon and lat should also not be processed unless required
    #xcoord = int(info[2])
    #ycoord = int(info[3])

    if int(info[3]) % dy_ts == 0 and int(info[2]) % dx_ts == 0 \
            and X_BND_PAD <= int(info[2]) < nx - X_BND_PAD \
            and Y_BND_PAD <= int(info[3]) < ny - Y_BND_PAD:
        name = hex(int('%s%s' % (info[2], info[3].zfill(4))))[2:].zfill(7).upper()
        llp.write('%10.4f %10.4f %s\n' \
                % (float(info[0]), float(info[1]), name))
        sc.append('%5d %5d %5d %s\n' \
                % (int(info[2]), int(info[3]), 1, name))

with open('%s.statcords' % (STATGRID_GEN), 'w') as scp:
    scp.write('%d\n' % (len(sc)))
    scp.write(''.join(sc))
llp.close()


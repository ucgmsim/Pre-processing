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

from shared import *
from params import *

verify_files([MODEL_COORDS])
verify_logfiles([STATGRID_GEN])
verify_strings([nx, ny, dx_ts, dy_ts, X_BND_PAD, Y_BND_PAD])

nx, ny, dx_ts, dy_ts = map(int,[nx, ny, dx_ts, dy_ts])

op = open(STATGRID_GEN, 'w')
mcp = open(MODEL_COORDS, 'r')
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
        op.write('%10.4f %10.4f %s\n' \
                % (float(info[0]), float(info[1]), \
                hex(int('%s%s' % (info[2], info[3])))[2:].zfill(7).upper()))

op.close()
mcp.close()


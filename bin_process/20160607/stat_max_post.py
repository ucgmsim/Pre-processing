#!/usr/bin/env python2

# h (0.0->1.0), s (0.0->1.0), v (0.0->255)
from colorsys import hsv_to_rgb
from shared import *
from params import *

round_int_str = lambda x: str(int(round(x)))

verify_files([MAX_STATFILE])

# find max for scaling factor
max_stat = 0.0
with open(MAX_STATFILE, 'r') as sp:
    for line in sp:
        if float(line.split()[2]) > max_stat:
            max_stat = float(line.split()[2])

if max_stat > 80.0:
    max_stat = 80.0

# modify lines before writing them back
old_lines = []
with open(MAX_STATFILE, 'r') as sp:
    old_lines = sp.readlines()

# Scaling Background:
# HSV: Saturation = 100%, Value/Lightness = 50%
# H: 1/3 (green minimum) -> 0 (red maximum)
# Motion: 0 (green minimum) -> max_stat (red maximum)
factor = 2.0/3.0 / max_stat
# 1/3 - (station max * factor) = desired hue

with open(MAX_STATFILE, 'w') as op:
    # GMT colour format: R/G/B@A (A: 0 opaque -> 100 transparent)
    for line in old_lines:
        if float(line.split()[2]) <= max_stat:
            h = 2.0/3.0 - (float(line.split()[2]) * factor)
            op.write('%s %s %s@20\n' % (line.rstrip(), str(max_stat), '/'.join(map(round_int_str, hsv_to_rgb(h, 1.0, 255.0)))))
        else:
            op.write('%s %s %s@20\n' % (line.rstrip(), str(max_stat), '255/0/255@20'))




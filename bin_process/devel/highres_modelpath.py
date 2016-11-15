#!/usr/bin/env python2

from params import *
from tools import *


with open(vel_mod_params, 'r') as mp:
    lines = mp.readlines()

corners = []
labels = ['c1=', 'c2=', 'c3=', 'c4=']

# find corners using labels
for label in labels:
    for line in lines:
        if label in line:
            corners.append(map(float, line.split()[1:3]))

# close the box by going back to origin
corners.append(corners[0])
# until each side has at least 100 points
while len(corners) < 4 * 100:
    # work backwards, insertions don't change later indexes
    for i in xrange(len(corners) - 1, 0, -1):
        val = ll_mid(corners[i][0], corners[i][1], \
                corners[i - 1][0], corners[i - 1][1])
        corners.insert(i, val)

with open('sim.modelpath_hr', 'w') as mp:
    for point in corners:
        mp.write('%s %s\n' % (point[0], point[1]))

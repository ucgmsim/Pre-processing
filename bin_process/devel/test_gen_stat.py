#!/usr/bin/env python2

from glob import glob
import os

from tools import ll_dist

files = glob('velBB/*.ver')
stats = []
for f in files:
    stat = f.split('/')[1].split('.')[0]
    stats.append(stat)

spaced = []
spaced.append([172.61990, -43.52930, 'CBGS'])
spaced.append([174.13840, -41.82740, 'WDFS'])
spaced.append([174.76818, -41.28405, 'WEL'])
spaced.append([173.05360, -42.61950, 'WTMC'])

with open('geonet_stations_20161119.ll', 'r') as gnf:
    for line in gnf:
        x, y, n = line.split()
        if n in stats:
            good = True
            for s in spaced:
                if ll_dist(float(x), float(y), s[0], s[1]) < 50:
                    good = False
                    break
            if good:
                spaced.append([float(x), float(y), n])

out = open('xy.ll', 'w')
for s in spaced:
    out.write('%s %s %s\n' % (s[0], s[1], s[2]))
out.close()

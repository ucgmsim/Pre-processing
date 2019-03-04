#!/usr/bin/env python2

from qcore import geo

nhm = 'NZ_FLTmodel_2010.txt'
with open(nhm, 'r') as nf:
    db = nf.readlines()
    dbi = 15
    dbl = len(db)
    while dbi < dbl:
        name = db[dbi].strip()
        n_pt = int(db[dbi + 11])
        pts = [list(map(float, ll.split())) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
        dists = [geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) \
                for i in range(n_pt - 1)]
        if min(dists) < 1:
            print(name)
            print('-' * 10)
            print(dists)
            print()
        dbi += 13 + n_pt

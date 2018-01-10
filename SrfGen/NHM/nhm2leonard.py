#!/usr/bin/env python2

import math

from qcore import geo

nhm = 'NZ_FLTmodel_2010.txt'

# Leonard 2014 Relations
def leonard(rake, A):
    # if dip slip else strike slip
    if round(rake % 360 / 90.) % 2:
        return 4.00 + math.log10(A)
    else:
        return 3.99 + math.log10(A)

# output header
print('Name,NHM_Mw,Leonard_Mw')

# loop through faults
with open(nhm, 'r') as nf:
    db = nf.readlines()
    dbi = 15
    dbl = len(db)
while dbi < dbl:
    name = db[dbi].strip()
    n_pt = int(db[dbi + 11])

    # basic properties
    Mw = float(db[dbi + 10].split()[0])
    dip = float(db[dbi + 3].split()[0])
    rake = float(db[dbi + 5])

    # find total area
    pts = [map(float, ll.split()) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
    dbottom = float(db[dbi + 6].split()[0])
    dtop = float(db[dbi + 7].split()[0])
    # Karim: add 3km to depth if bottom >= 12km
    extend = dbottom >= 12
    # fault width (along dip)
    fwid = (dbottom - dtop + 3 * extend) / math.sin(math.radians(dip))
    flen = sum([geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) \
            for i in xrange(n_pt - 1)])
    Ml = leonard(rake, fwid * flen)

    # output
    print('%s,%s,%s' % (name, Mw, Ml))

    # move to next definition
    dbi += 13 + n_pt

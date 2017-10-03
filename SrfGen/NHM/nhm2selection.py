#!/usr/bin/env python2

import math
import os
import sys
from tempfile import mkstemp

from subprocess import call

import geo
import gmt

NHM_FILE = 'NZ_FLTmodel_2010.txt'
NHM_START = 15
N_HYPO = '20k'
N_SLIP = 3
# 'res/ni_fault_range.txt', 'res/si_fault_range.txt'
CLIP_POLYGON = 'res/ni_fault_range.txt'
OUT_FILE = 'nhm_selection_geo.txt'

out_dir = os.path.dirname(OUT_FILE)
if out_dir == '':
    out_dir = '.'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
out = open(OUT_FILE, 'w')

# load db
dbi = NHM_START
with open(NHM_FILE, 'r') as dbr:
    db = map(str.strip, dbr.readlines())
dbl = len(db)

# find criteria matching faults
while dbi < dbl:
    ###
    ### parameters from db
    ###
    name = db[dbi]
    dip = float(db[dbi + 3].split()[0])
    dip_dir = float(db[dbi + 4])
    dbottom = float(db[dbi + 6].split()[0])
    dtop = float(db[dbi + 7].split()[0])
    n_pt = int(db[dbi + 11])
    pts = [map(float, ll.split()) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
    # Karim: add 3km to depth if bottom >= 12km
    extend = dbottom >= 12
    # fault plane width (along dip)
    fwid = (dbottom - dtop + 3 * extend) \
            / math.sin(math.radians(dip))
    # projected (surface) width
    pwid = fwid * math.cos(math.radians(dip))
    # rest of surface trace
    pts.extend([geo.ll_shift(ll[1], ll[0], pwid, dip_dir)[::-1] \
            for ll in pts[::-1]])

    ###
    ### find if part of fault in south island
    ###
    # store trace in file
    trace_file = mkstemp()
    try:
        with open(trace_file[1], 'w') as tf:
            tf.write('\n'.join([' '.join(map(str, ll)) for ll in pts]))
        if len(gmt.truncate(trace_file[1], clip = CLIP_POLYGON)):
            print('%s %s %s' % (name, N_HYPO, N_SLIP))
            out.write('%s %s %s\n' % (name, N_HYPO, N_SLIP))
    except Exception as e:
        print(e)
        sys.exit(1)
    finally:
        os.remove(trace_file[1])

    # work on next fault
    dbi += 13 + n_pt

out.close()

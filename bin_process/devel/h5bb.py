#!/usr/bin/env python2

from glob import glob
from multiprocessing import Pool
from os.path import basename
import sys

import h5py as h5
import numpy as np

from params_base_bb import *

procs = 64
if len(sys.argv) > 1:
    procs = int(sys.argv[1])

COMP_EXTS = {'090':0, '000':1, 'ver':2}
file_list = glob('%s/*_???????.???.bb' % (hf_accdir))
num_files = len(file_list)
if num_files < procs:
    procs = num_files
proc_lists = [file_list[g::procs] for g in xrange(procs)]

def write2hdf(proc_list):
    h5p = h5.File('virtual.hdf5', 'a')
    for ii, bbf in enumerate(proc_list):
        print '%d of %d' % (ii + 1, num_files)
        # read binary values
        bb_vel = np.fromfile(bbf, '>f')

        # read from intermediate binary
        f_details = basename(bbf).split('.')
        comp = COMP_EXTS[f_details[1]]
        xy = str(int(f_details[0].split('_')[-1], 16)).zfill(8)
        x = xy[:4]
        y = xy[4:]
        # write back to HDF5 file
        data = h5p['%s%s/VEL' % (x, y)][:, comp] = bb_vel

    # close HDF5
    h5p.close()

p = Pool(procs)
p.map(write2hdf, proc_lists)

#!/usr/bin/env python2

from glob import glob
from os.path import basename

import h5py as h5
import numpy as np

COMP_EXTS = {'090':0, '000':1, 'ver':2}
h5p = h5.File('virtual.hdf5', 'a')

file_list = glob('RunFolder/acc/*.*')
num_files = len(file_list)

for ii, bbf in enumerate(file_list):
    print '%d of %d' % (ii + 1, num_files)
    # read binary values
    bb_vel = np.fromfile(hff, '>f')

    # read from HDF5
    f_details = basename(bbf).split('.')
    comp = COMP_EXTS[f_details[1]]
    xy = str(int(f_details[0].split('_')[-1], 16)).zfill(8)
    x = xy[:4]
    y = xy[4:]
    # write back to HDF5 file
    h5p['%s%s/VEL' % (x, y)][..., comp] = bb_vel

# close HDF5
h5p.close()


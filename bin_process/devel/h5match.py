#!/usr/bin/env python2

from glob import glob
from multiprocessing import Pool
from os import remove
from os.path import basename
# math functions faster than numpy for simple data
from math import ceil, log, pi, tan, sqrt

import h5py as h5
import numpy as np
rfft = np.fft.rfft
irfft = np.fft.irfft
# sosfilt new in scipy 0.16
# sosfiltfilt new in scipy 0.19
from scipy.signal import butter, sosfilt
import matplotlib.pyplot as plt

from siteamp_models import cb08_amp

# TODO: read params from params file
procs = 4
nt = 20000
dt = 0.005
ft_len = int(2 ** ceil(log(nt)/log(2)))
nyq = 1.0 / (2.0 * dt)
COMP_EXTS = {'090':0, '000':1, 'ver':2}
# frequencies in the fourier domain
#freqs = 200.0 * (np.arange(0, 32768/2) / 32768.0)
file_list = glob('RunFolder/acc/*.*')
proc_lists = []
for g in procs:
    proc_lists.append([file_list[g::procs]])

# butterworth filter
# bandpass not necessary as sampling frequency too low
def bwfilter(data, freq, band):
    """
    data: np.array to filter
    freq: cutoff frequency
    band: 'highpass' or 'lowpass'
    """
    return sosfilt(butter(4, freq / nyq, btype = band, output = 'sos'), data)

def run_match(proc_list):
    num_files = len(proc_list)
    h5p = h5.File('virtual.hdf5', 'r')
    for ii, hff in enumerate(proc_list):
        print '%d of %d' % (ii + 1, num_files)
        ###
        # PROCESS HF

        # read binary values
        v = np.fromfile(hff, '>f')

        # if using ascii input
        #v = np.empty(nt)
        #with open(hff, 'rb') as f:
        #    l1 = f.readline()
        #    l2 = f.readline()
        #    for i, s in enumerate(map(float, \
        #                ' '.join(map(str.rstrip, f.readlines())).split())):
        #        v[i] = s

        # peak ground acceleration
        pga = np.max(np.abs(v)) / 981.0
        # taper ntap values on the right using the hanning method
        ntap = nt * 0.05
        v *= dt
        v[nt - ntap:] *= np.hanning(ntap * 2 + 1)[ntap + 1:]

        # extend with blanks for fft
        v.resize(ft_len)
        ft = rfft(v)
        # amplification factors are applied on fft values
        ampf = cb08_amp(0.005, ft_len, 865.0, 250.0, 850.0, pga)
        ft[:-1] *= ampf
        s = irfft(ft)[:nt] * (nt * dt * 2)

        # filter unwanted frequencies
        hf_acc = bwfilter(s, 1.0, 'highpass')

        ###
        # PROCESS LF

        # read from HDF5
        f_details = basename(hff).split('.')
        comp = COMP_EXTS[f_details[1]]
        xy = str(int(f_details[0].split('_')[-1], 16)).zfill(8)
        x = xy[:4]
        y = xy[4:]
        try:
            lf_vel = h5p['%s%s/VEL' % (x, y)][..., comp]
        except KeyError:
            # this station doesn't exist in LF sim
            remove(hff)
            print('WARNING: station at (%d, %d) not found for LF sim.' % (x, y))
            continue
        lf_acc = np.diff(np.hstack(([0], lf_vel))) * (1.0 / dt)
        # fft
        lf_acc.resize(ft_len)
        lf_ft = rfft(lf_acc)
        # amp
        lf_ft[:-1] *= ampf
        # ifft
        lf_acc = irfft(lf_ft)[nt]
        # filter
        lf_acc = bwfilter(s, 1.0, 'lowpass')

        ###
        # PROCESS BB

        # add values together and integrate ACC to VEL
        bb_vel = np.cumsum((hf_acc + lf_acc) * dt)

        # write back to file
        bb_vel.tofile(hff)

    # close HDF5
    h5p.close()


# run each process with a different set of files to work on
p = Pool(procs)
p.map(run_match, proc_lists)


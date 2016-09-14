"""
Shared functions to work on time-series.

@author Viktor Polak
@date 13/09/2016
"""

from math import ceil, log

import numpy as np
rfft = np.fft.rfft
irfft = np.fft.irfft

from siteamp_models import cb08_amp

def get_ft_len(nt):
    """
    Length the fourier transform should be
    given timeseries length nt.
    """
    return int(2 ** ceil(log(nt) / log(2)))

def ampdeamp(timeseries, ampf, amp = True):
    """
    Amplify or Deamplify timeseries.
    """
    nt = len(timeseries)

    # length the fourier transform should be
    ft_len = get_ft_len(nt)

    # taper 5% on the right using the hanning method
    ntap = int(nt * 0.05)
    timeseries[nt - ntap:] *= np.hanning(ntap * 2 + 1)[ntap + 1:]

    # extend array, fft
    timeseries = np.resize(timeseries, ft_len)
    timeseries[nt:] = 0
    fourier = rfft(timeseries)

    # ampf modified for de-amplification
    if not amp:
        ampf = 1 / ampf
    # last value of fft is some identity value
    fourier[:-1] *= ampf

    return irfft(fourier)[:nt]

def read_ascii(filepath):
    """
    Read timeseries data from standard ascii file to numpy array.
    """
    with open(filepath, 'r') as ts:
        info1 = ts.readline().split()
        info2 = ts.readline().split()

        nt = int(info2[0])
        vals = np.empty(nt)
        for i, a_val in enumerate(map(float, \
                ' '.join(map(str.rstrip, ts.readlines())).split())):
            vals[i] = a_val

    # check if header length correct
    # otherwise end of array may be garbage
    assert(i == nt - 1)

    return vals

def vel2acc(timeseries, dt):
    """
    Differentiate following Rob Graves' code logic.
    """
    return np.diff(np.hstack(([0], timeseries)) * (1.0 / dt))

def acc2vel(timeseries, dt):
    """
    Integrates following Rob Graves' code logic (simple).
    """
    return np.cumsum(timeseries) * dt


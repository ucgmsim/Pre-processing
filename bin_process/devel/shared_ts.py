"""
Shared functions to work on time-series.

@author Viktor Polak
@date 13/09/2016
"""

from math import ceil, log, pi

# sosfilt new in scipy 0.16
# sosfiltfilt new in scipy 0.18
try:
    from scipy.signal import butter
    # this is a local import
    from sosfiltfilt import sosfiltfilt
except ImportError:
    print('SciPy not installed. Certain functions will fail.')
import numpy as np
rfft = np.fft.rfft
irfft = np.fft.irfft

from params import *

# butterworth filter
# bandpass not necessary as sampling frequency too low
def bwfilter(data, dt, freq, band):
    """
    data: np.array to filter
    freq: cutoff frequency
    band: 'highpass' or 'lowpass'
    """
    # power spectrum based LF/HF filter (shift cutoff)
    # readable code commented, fast code uncommented
    #order = 4
    #x = 1.0 / (2.0 * order)
    #if band == 'lowpass':
    #    x *= -1
    #freq *= exp(x * log(sqrt(2.0) - 1.0))
    nyq = 1.0 / (2.0 * dt)
    if match_powersb:
        if band == 'highpass':
            freq *= 0.8956803352330285
        else:
            freq *= 1.1164697500474103
    return sosfiltfilt( \
            butter(4, freq / nyq, btype = band, output = 'sos'), \
            data, padtype = None)

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

def transf(vs_soil, rho_soil, damp_soil, height_soil, \
        vs_rock, rho_rock, damp_rock, nt, dt):
    """
    Used in deconvolution. Made by Chris de la Torre.
    vs = shear wave velocity (upper soil or rock)
    rho = density
    damp = damping ratio
    height_soil = height of soil above rock
    nt = number of timesteps
    dt = delta time in timestep (seconds)
    """
    ft_len = get_ft_len(nt)
    # TODO: before it was ft_len / 2 + 1 but this may be an error
    # the last value isn't an ft value
    ft_freq = np.arange(0, ft_len / 2) * (1 / (ft_len * dt))

    omega = 2.0 * pi * ft_freq
    Gs = rho_soil * vs_soil ** 2.0
    Gr = rho_rock * vs_rock ** 2.0

    kS = omega / (vs_soil * (1.0 + 1j * damp_soil))
    kR = omega / (vs_rock * (1.0 + 1j * damp_rock))

    alpha = Gs * kS / (Gr * kR)

    H = 2.0 / ((1.0 + alpha) * np.exp(1j * jS * hS) + (1.0 - alpha) \
            * np.exp(-1j * kS * hS))
    H[0] = 1

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


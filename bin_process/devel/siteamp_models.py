"""
Acceleration amplification models.

@date 24 June 2016
@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz

Implemented Models
==============================
cb08_amp: Based on Campbell and Bozorgnia 2008 - added 24 June 2016

Usage
==============================
from siteamp_models import cb08_amp
cb08_amp(variables, ...)
"""

# math functions faster than numpy for non-vector data
from math import ceil, exp, log

import numpy as np
# divide by zero = inf
np.seterr(divide='ignore')


# Campbell Bozorgnia 2008
# define constant values outside function to prevent re-calculation
# TODO: prefix variable names with 'cb08_' with multiple implementations
n_per = 22
scon_c = 1.88
scon_n = 1.18
per = np.array([0.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, \
        0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, \
        2.00, 3.00, 4.00, 5.00, 7.50, 10.0])
m_freq = 1.0 / per
f1_src = np.hstack(([1000.0], m_freq[1:]))
c10 = np.array([1.058, 1.058, 1.102, 1.174, 1.272, 1.438, 1.604, 1.928, \
        2.194, 2.351, 2.460, 2.587, 2.544, 2.133, 1.571, 0.406, \
        -0.456, -0.82, -0.82, -0.82, -0.82, -0.82])
c11 = np.array([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, \
        0.04, 0.04, 0.04, 0.04, 0.04, 0.077, 0.15, 0.253, \
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
k1 = np.array([865.0, 865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, 878.0, \
        748.0, 654.0, 587.0, 503.0, 457.0, 410.0, 400.0, 400.0, \
        400.0, 400.0, 400.0, 400.0, 400.0, 400.0])
k2 = np.array([-1.186, -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, \
        -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025,
        -0.299, 0.0, 0.0, 0.0, 0.0, 0.0])
k3 = np.array([1.839, 1.839, 1.840, 1.841, 1.843, 1.845, 1.847, 1.852, \
        1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, \
        2.019, 2.110, 2.200, 2.291, 2.517, 2.744])
flowcap = 0.0
fmin = 0.05
fmidbot = 0.1
fmid = fmidbot
fhigh = fmid
fhightop = 20.0
fmax = 50.0
# f_site function domains
fs_low = lambda T, vs30, a1100 : c10[T] * log(vs30 / k1[T]) + \
        k2[T] * log((a1100 + scon_c * exp(scon_n * log(vs30 / k1[T]))) / (a1100 + scon_c))
fs_mid = lambda T, vs30, a1100 = None : (c10[T] + k2[T] * scon_n) * log(vs30 / k1[T])
fs_high = lambda T, vs30 = None, a1100 = None : (c10[T] + k2[T] * scon_n) * log(1100.0 / k1[T])
fs_auto = lambda T, vs30 : fs_low if vs30 < k1[T] else fs_mid if vs30 < 1100.0 else fs_high
fs1100 = fs_high(0)

def cb08_amp(dt, n, vref, vsite, vpga, pga):
    # default amplification is 1.0 (keeping values the same)
    ampf = np.ones(n / 2, np.float)

    fs_vpga = fs_auto(0, vpga)(0, vpga, pga)
    a1100 = pga * exp(fs1100 - fs_vpga)

    # calculate factor for each period
    it = (exp(fs_auto(T, vsite)(T, vsite, a1100) - fs_auto(T, vref)(T, vref, a1100)) \
            for T in xrange(n_per))
    ampf0 = np.fromiter(it, np.float, count = n_per)

    try:
        # T is the first occurance of a value <= flowcap
        # throws IndexError if no results (the second [0])
        T = np.nonzero((m_freq[:-1] <= flowcap))[0][0]
        # T cannot be the last value because of the following logic
        ampf0[T + 1:] = ampf0[T]
    except IndexError:
        pass
    # frequencies of fourier transform
    ftfreq = np.arange(1, n / 2) * (1.0 / (n * dt))
    # TODO: vectorise to improve speed
    #a0 = np.repeat(ampf0[-1], len(ftfreq))
    #f0 = np.repeat(f1_src[-1], len(ftfreq))
    j = n_per - 1
    f0 = f1_src[j]
    a0 = ampf0[j]
    f1 = f0
    a1 = a0
    dadf = 0.0
    for i, ftf in enumerate(ftfreq):
        if ftf > f1:
            f0 = f1
            a0 = a1
            if j - 1:
                j -= 1
                a1 = ampf0[j]
                f1 = f1_src[j]
                dadf = (a1 - a0) / log(f1 / f0)
            else:
                dadf = 0.0
    #dadf = (a1 - a0) / np.log(f1 / f0) * (f1 == f0)
        ampv = a0 + dadf * log(ftf / f0)

        if ftf < fmin:
            continue
        if ftf < fmidbot:
            ampf[i + 1] = 1.0 + log(ftf / fmin) * (ampv - 1.0) / log(fmidbot / fmin)
        elif ftf < fhightop:
            ampf[i + 1] = ampv
        elif ftf < fmax:
            ampf[i + 1] = ampv + log(ftf / fhightop) * (1.0 - ampv) / log(fmax / fhightop)

    return ampf


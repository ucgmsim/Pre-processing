#!/usr/bin/env python2

# math functions faster than numpy for simple data
from math import ceil, exp, log

import numpy as n

# Campbell Bozorgnia 2008
# define constant values outside function to prevent re-calculation
n_per = 22
scon_c = 1.88
scon_n = 1.18
per = n.array([0.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, \
        0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, \
        2.00, 3.00, 4.00, 5.00, 7.50, 10.0])
freq = 1.0 / per
c10 = n.array([1.058, 1.058, 1.102, 1.174, 1.272, 1.438, 1.604, 1.928, \
        2.194, 2.351, 2.460, 2.587, 2.544, 2.133, 1.571, 0.406, \
        -0.456, -0.82, -0.82, -0.82, -0.82, -0.82])
c11 = n.array([0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, \
        0.04, 0.04, 0.04, 0.04, 0.04, 0.077, 0.15, 0.253, \
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03])
k1 = n.array([865.0, 865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, 878.0, \
        748.0, 654.0, 587.0, 503.0, 457.0, 410.0, 400.0, 400.0, \
        400.0, 400.0, 400.0, 400.0, 400.0, 400.0])
k2 = n.array([-1.186, -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, \
        -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025,
        -0.299, 0.0, 0.0, 0.0, 0.0, 0.0])
k3 = n.array([1.839, 1.839, 1.840, 1.841, 1.843, 1.845, 1.847, 1.852, \
        1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, \
        2.019, 2.110, 2.200, 2.291, 2.517, 2.744])
# f_site function domains
fs_low = lambda T, vs30, a1100 : c10[T] * log(vs30 / k1[T]) + \
        k2[T] * log(a1100 + scon_c * exp(scon_n * log(vs30 / k1[T])) / (a1100 + scon_c))
fs_mid = lambda T, vs30, a1100 = None : (c10[T] + k2[T] * scon_n) * log(vs30 / k1[T])
fs_high = lambda T, vs30 = None, a1100 = None : (c10[T] + k2[T] * scon_n) * log(1100.0 / k1[T])
fs_auto = lambda vs30 : fs_low if vs30 < k1[0] else fs_mid if vs30 < 1100 else fs_high
fs1100 = fs_high(0)
def cb2008_ampf(ampf, dt, n, vref, vsite, vpga, pga, fmin, fmidbot, fmid, fhigh, fhightop, fmax, flowcap):
    fs_vpga = fs_auto(vpga)(0, vpga, pga)
    a1100 = pga * exp(fs1100 - fs_vpga)
    it = (exp(fs_auto(vsite)(T, vsite, a1100) - fs_auto(vref)(T, vref, a1100)) \
            for T in xrange(n_per))
    ampf0 = n.fromiter(it, np.float, count = n_per)
    try:
        # T is the first occurance of a value <= flowcap
        # throws IndexError if no results (the second [0])
        T = n.nonzero((freq[:-1] <= flowcap))[0][0]
        # T cannot be the last value because of the following logic
        ampf0[T + 1:] = ampf0[T]
    except IndexError:
        pass


nt = 20000
dt = 0.005

getnt_p2 = lambda nt : int(2 ** ceil(log(nt)/log(2)))

f = open('test.ver', 'r')
f.readline()
f.readline()

v = n.zeros(nt)
for i, s in enumerate(map(float, \
            ' '.join(map(str.rstrip, f.readlines())).split())):
    v[i] = s

ntap = nt * 0.05
v *= dt
v[nt - ntap:] *= n.hanning(ntap * 2 + 1)[ntap + 1:]

l = getnt_p2(nt)
print len(v)
print len(n.fft.rfft(v))
v.resize(l)
ft = n.fft.rfft(v)
ft2 = n.dstack((ft.real, ft.imag)).flatten()
ft2[1] = ft2[-2]
print ft2[:l]

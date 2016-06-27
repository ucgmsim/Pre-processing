#!/usr/bin/env python2

# math functions faster than numpy for simple data
from math import ceil, exp, log

import numpy as np

from siteamp_models import cb08_amp

nt = 20000
dt = 0.005

getnt_p2 = lambda nt : int(2 ** ceil(log(nt)/log(2)))

f = open('test.ver', 'r')
f.readline()
f.readline()

v = np.zeros(nt)
for i, s in enumerate(map(float, \
            ' '.join(map(str.rstrip, f.readlines())).split())):
    v[i] = s
pga = np.max(np.abs(v))/981.0

ntap = nt * 0.05
v *= dt
v[nt - ntap:] *= np.hanning(ntap * 2 + 1)[ntap + 1:]

l = getnt_p2(nt)

v.resize(l)
ft = np.fft.rfft(v)
#print v[4800:4810]
#print np.fft.irfft(ft)[4800:4810]
#ft2 = np.dstack((ft.real, ft.imag)).flatten()
#ft2[1] = ft2[-2]
#ft2 = ft2[:l]

ampf = cb08_amp(0.005, l, 865.0, 250.0, 850.0, pga, 0.5, 0.5, 1.0, 3.333, 10.0, 15.0, 0.0)

rft = ft.real
rft[:-1] *= ampf
ift = ft.imag
ift[:-1] *= ampf
#ift[0] = ft.real[-1]
#ift *= ampf

# ^^ working
fta = (rft + 1j * ift)
#fta = np.hstack((fta,ft[-1] * ampf[-2]))
real = fta[0].real
imag = fta[0].imag
#fta[0] = real + imag + 1j * (real - imag)
s = np.fft.irfft(fta) * 200

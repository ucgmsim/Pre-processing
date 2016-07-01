#!/usr/bin/env python2
"""
Displays summary of accelerograms for human evaluation.

@author Viktor Polak
@contact viktor.polak@canterbury.ac.nz
@date 14 June 2016

ISSUES: 
"""

from glob import glob
from math import sin, cos, radians

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np

# line length (values per line)
VPL = 10
# lines in header
SIZE_HEAD = 26

# converts bearing format (N|S)0-90(E|W) to degrees clockwise from north
# sample input: "S28E"
def nsew2deg(fstring):
    if fstring[0] == 'N':
        if fstring[-1] == 'E':
            return int(fstring[1:-1])
        else:
            return 360 - int(fstring[1:-1])
    if fstring[0] == 'S':
        if fstring[-1] == 'W':
            return 180 + int(fstring[1:-1])
        else:
            return 180 - int(fstring[1:-1])
    

file_pattern = 'Vol1/data/*.V1A'
inputs = glob(file_pattern)
# testing - only work on first file
#inputs = inputs[:1]
pdf = PdfPages('plots.pdf')

for vfile in inputs:
    with open(vfile, 'r') as vp:
        data = map(str.rstrip, vp.readlines())

    # following is assumed to be the same amongst components
    stat = data[1].split()[1]       # station name (short form)
    stat_h = data[2]                # human readable station name

    # these components can be dealt with if different amongst components
    nts = []    # number of time steps
    dts = []    # time step
    vls = []    # number of text lines containing data
    cns = []    # component names

    # read components seperately
    # uncleansed files can have differing params such as nt
    # jump is the start line position of the component in the file
    jump = 0
    while jump >= 0:
        nts.append(int(data[jump + 19].split()[0]))   # number of timesteps
        dts.append(float(data[jump + 22].split()[4])) # timestep
        vl = nts[-1] // VPL                     # lines containing values
        if nts[-1] % VPL != 0:
            vl += 1
        vls.append(vl)
        cns.append(data[jump + 12].split()[1])
        jump = jump + SIZE_HEAD + vl
        if jump >= len(data):
            jump = -1

    fig = plt.figure()
    fig.suptitle('[%s] %s' % (stat, stat_h), fontsize=18, fontweight='bold')
    for i in xrange(len(cns)):
        # x = time
        # numpy thinks that np.arrange(0, 64.805, 0.005) includes 64.805 (rounding?)
        # therefore takes 0.1 times the timestep away to prevent this behaviour
        x = np.arange(0, nts[i] * dts[i] - 0.1 * dts[i], dts[i])
        # content
        ax = plt.subplot(len(cns) * 100 + 10 + i + 1)
        ax.set_title(cns[i])
        dstart = SIZE_HEAD * i + sum(vls[:i]) + 26
        y = np.array(map(float, ' '.join(data[dstart:dstart + vls[i]]).split()))
        ax.plot(x, np.array(map(float, ' '.join( \
                data[dstart:dstart + vls[i]]).split())), linewidth=0.2)
        plt.xlim(0, 100)
        # using subplots will cause too many ticks
        plt.locator_params(axis = 'y', nbins = 4)
    # prevest titles and axis lables from overlapping
    fig.tight_layout()
    # leave room for title
    plt.subplots_adjust(top = 0.85)
    pdf.savefig()
    plt.close()
    
    # only used if correcting

    #if (c1 + 90) % 360 == c2:
    #    # to be 000 component is first
    #    theta = radians(c1)
    #    flipc12 = False
    #else:
    #    # to be 000 component is second
    #    theta = radians(c2)
    #    flipc12 = True

    
    # rotates components, to invert vertical component, change 3,3 from 1 to -1
    #rot_m = np.array([[cos(theta),  -sin(theta), 0], \
    #                  [-sin(theta), -cos(theta), 0], \
    #                  [0,           0,           1]])

    #comps = np.dot(np.dstack(( \
    #    np.array(map(float, ' '.join(data[clen * 0 + 26:clen * 0 + 26 + vlines]).split())), \
    #    np.array(map(float, ' '.join(data[clen * 1 + 26:clen * 1 + 26 + vlines]).split())), \
    #    np.array(map(float, ' '.join(data[clen * 2 + 26:clen * 2 + 26 + vlines]).split())) \
    #    )), rot_m)[0]

pdf.close()

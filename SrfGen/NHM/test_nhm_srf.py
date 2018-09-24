#!/usr/bin/env python

from argparse import ArgumentParser
from glob import glob
import math
from multiprocessing import Pool
import os
from shutil import copy, rmtree
from subprocess import call
from tempfile import mkdtemp

from h5py import File as h5open
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np

from createSRF import leonard
from qcore.geo import ll_shift, ll_bearing
from qcore import gmt

def test_mw_vs_area(srf_infos):
    # load mw and areas from srf infos
    ds_mw = []
    ds_area = []
    ss_mw = []
    ss_area = []
    for i in srf_infos:
        with h5open(i, 'r') as h:
            area = sum(h.attrs['length'] * h.attrs['width'])
            # if dip slip else strike slip
            if round(h.attrs['rake'][0][0] % 360 / 90.) % 2:
                ds_mw.append(h.attrs['mag'])
                ds_area.append(area)
            else:
                ss_mw.append(h.attrs['mag'])
                ss_area.append(area)

    # plot
    fig = plt.figure(figsize = (8, 5), dpi = 150)
    min_y = 10
    max_y = 0
    if len(ds_mw):
        min_y = min(min(ds_mw), min_y)
        max_y = max(max(ds_mw), max_y)
    if len(ss_mw):
        min_y = min(min(ss_mw), min_y)
        max_y = max(max(ss_mw), max_y)
    y = np.linspace(min_y, max_y, 100)

    if len(ds_mw):
        plt.plot(10**(y-4), y, label='Leonard (dip slip)')
        plt.plot(ds_area, ds_mw, label='SRF (dip slip)', marker='x', linestyle='None')
    if len(ss_mw):
        plt.plot(10**(y-3.99), y, label='Leonard (strike slip)')
        plt.plot(ss_area, ss_mw, label='SRF (strike slip)', marker='x', linestyle='None')

    plt.legend(loc='best')
    plt.gca().set_xscale('log')
    plt.title('Area vs Magnitude')
    plt.ylabel('Magnitude')
    plt.xlabel('Area (sq. km)')
    plt.show()
    plt.close()

    # part 2 : mw cs vs mw leonard
    fig = plt.figure(figsize = (8, 5), dpi = 150)
    plt.plot(y, y, label='Target')
    if len(ds_mw):
        plt.plot(ds_mw, 4 + np.log10(ds_area), label='Actual (dip slip)', marker='x', linestyle='None')
    if len(ss_mw):
        plt.plot(ss_mw, 3.99 + np.log10(ss_area), label='Actual (strike slip)', marker='x', linestyle='None')
    plt.legend(loc='best')
    plt.title('Mw SRF vs Mw Leonard')
    plt.xlabel('Mw SRF')
    plt.ylabel('Mw Leonard')
    plt.show()
    plt.close()

def test_mw_vs_nrup(srf_dirs):
    nrup = []
    mw = []
    for d in srf_dirs:
        infos = glob(os.path.join(d, '*.info'))
        nrup.append(len(infos))
        with h5open(infos[0], 'r') as h:
            mw.append(h.attrs['mag'])

    # plot
    x = []
    x.append(min(mw))
    x.append(max(mw))
    if x[0] < 6 < x[-1]:
        x.insert(1, 6)
    if x[0] < 8 < x[-1]:
        x.insert(-1, 8)

    fig = plt.figure(figsize = (8, 5), dpi = 150)
    plt.plot(x, np.maximum(np.minimum((np.array(x) * 20 - 110), 50), 10), label = "Target")
    plt.plot(mw, nrup, label='Actual', marker='x', linestyle='None')

    plt.legend(loc='best')
    plt.title('Magnitude vs Number of Realizations')
    plt.ylabel('Number of Realizations')
    plt.xlabel('Magnitude')
    plt.show()
    plt.close()

def test_seismogenic_depth(srf_dirs, nhm_file):
    srf_depths = []
    fault_names = []
    nhm_depths = []
    for d in srf_dirs:
        info = glob(os.path.join(d, '*.info'))[0]
        with h5open(info, 'r') as i:
            srf_depths.append(np.average(i.attrs['dbottom']))
        fault_names.append(os.path.basename(os.path.dirname(d)))
    # grab depths from nhm
    nhm_depths = np.zeros(len(fault_names))
    with open(nhm_file, 'r') as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
        n_i = 0
        while min(nhm_depths) <= 0:
            name = nhm[n_i].strip()
            if name in fault_names:
                nhm_depths[fault_names.index(name)] = float(nhm[n_i + 6].split()[0])
            elif name == '':
                print('Could not find all SRFs in NHM.')
                return
            n_i += 13 + int(nhm[n_i + 11])

    # plot
    fig = plt.figure(figsize = (8, 5), dpi = 150)
    max_x = max(max(srf_depths), max(nhm_depths))
    plt.plot([0, max_x], [0, max_x], label="NHM")
    plt.plot([0, max_x], [3, max_x + 3], label="NHM + 3")
    plt.plot(nhm_depths, np.array(srf_depths), label='SRF Depths', marker='x', linestyle='None')

    plt.legend(loc='best')
    plt.title('Seismogenic Depth from NHM and SRF')
    plt.ylabel('SRF Depth (km)')
    plt.xlabel('NHM Depth (km)')
    plt.show()
    plt.close()

def nhm2corners(nhm_file, names):
    corners_sub = []
    corners_shal = []
    ftypes = {}
    with open(nhm_file, 'r') as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
    n_i = 0
    while n_i < len(nhm):
        n = nhm[n_i].strip()
        n_pt = int(nhm[n_i + 11])
        ftype = 2 + int(not nhm[n_i + 1].startswith('SUBDUCTION'))
        dbottom = float(nhm[n_i + 6].split()[0])
        dtop = float(nhm[n_i + 7].split()[0])
        if n not in names:
            pwid = (dbottom - dtop + 3 * (dbottom >= 12)) \
                   / math.tan(math.radians(float(nhm[n_i + 3].split()[0])))
            ftype -= 2
            corners = np.zeros((n_pt - 1, 4, 2))
            trace = nhm[n_i + 12:n_i + 12 + n_pt]
            trace = [map(float, pair) for pair in map(str.split, trace)]
            dip_dir = float(nhm[n_i + 4])
            for i in range(len(corners)):
                corners[i, :2] = trace[i:i+2]
                corners[i, 2] = ll_shift(corners[i, 1, 1], corners[i, 1, 0], \
                                         pwid, dip_dir)[::-1]
                corners[i, 3] = ll_shift(corners[i, 0, 1], corners[i, 0, 0], \
                                         pwid, dip_dir)[::-1]
            if ftype % 2:
                corners_shal.append(corners)
            else:
                corners_sub.append(corners)
        ftypes[n] = ftype
        n_i += 13 + n_pt

    return ftypes, corners_shal, corners_sub

def test_spacial_srf(srf_dirs, nhm_file):
    tmp = mkdtemp()
    cpt = os.path.join(tmp, 'types.cpt')
    cpt_labels = {'0':';NHM Subduction\n', \
                  '1':';NHM Shallow\n', \
                  '2':';SRF Subduction\n', \
                  '3':';SRF Shallow\n'}
    gmt.makecpt('rainbow', cpt, 0, 4, inc=1, continuous=False)
    cpt_colours = []
    with open(cpt, 'r') as c:
        cpt_lines = c.readlines()
    for i in range(len(cpt_lines)):
        start = cpt_lines[i].split()[0]
        if start in cpt_labels:
            cpt_colours.append(cpt_lines[i].split()[1])
            cpt_lines[i] = '%s%s' % (cpt_lines[i].strip(), cpt_labels[start])
    with open(cpt, 'w') as c:
        c.write(''.join(cpt_lines))

    p = gmt.GMTPlot(os.path.join(tmp, 'srfs.ps'))
    p.spacial('M', region=gmt.nz_region, sizing=10, x_shift=1, y_shift=1)
    p.basemap(topo=None, road=None, highway=None, land='lightgray', water='white', res='f')
    p.ticks(major='2d')
    p.leave()
    copy(p.psf.name, os.path.join(tmp, 'srfs_ex.ps'))
    q = gmt.GMTPlot(os.path.join(tmp, 'srfs_ex.ps'), append=True, reset=False)
    p.enter()

    title = 'SRF and rest of NHM faults'
    zmax = 'NaN'
    length = 8 * 0.3
    for x in p, q:
        x.text(sum(gmt.nz_region[:2]) / 2.0, gmt.nz_region[3], \
               title, size='26p', dy=0.3)
        x.cpt_scale('C', 'R', cpt, pos='rel_out', horiz=False, arrow_f=False, \
                length=length, thickness=0.3, dx=0.3, categorical=True, gap='0.3', \
                zmax=zmax)
        title = 'Rest of NHM faults'
        zmax = 2
        length = 4 * 0.3

    def corners2gmt(corners):
        return '\n>\n'.join(['\n'.join([' '.join(map(str, ll)) for ll in plane]) for plane in corners])

    # nhm faults
    faults = map(os.path.basename, map(os.path.dirname, srf_dirs))
    ftypes, nhm_cnrs_shal, nhm_cnrs_sub = nhm2corners(nhm_file, faults)
    for x in p, q:
        x.path('\n>\n'.join(map(corners2gmt, nhm_cnrs_shal)), is_file=False, \
               colour='40/40/40', split='-', fill='%s@80' % (cpt_colours[1]), \
               close=True, width='0.2p')
        x.path('\n>\n'.join(map(corners2gmt, nhm_cnrs_sub)), is_file=False, \
               colour='40/40/40', split='-', fill='%s@80' % (cpt_colours[0]), \
               close=True, width='0.2p')
        x.path('\n>\n'.join([corners2gmt(c[:, :2]) for c in nhm_cnrs_shal]), \
               is_file=False, colour='40/40/40')
        x.path('\n>\n'.join([corners2gmt(c[:, :2]) for c in nhm_cnrs_sub]), \
               is_file=False, colour='40/40/40')

    # srf files
    for i, d in enumerate(srf_dirs):
        name = faults[i]
        f = glob(os.path.join(d, "*.info"))[0]
        with h5open(f, 'r') as h:
            corners = h.attrs['corners']
        # slower alternative if corners are bad
        #p.fault(glob(os.path.join(d, "*.srf"))[0], is_srf=True)
        p.path(corners2gmt(corners), is_file=False, close=True, split='-', fill='%s@80' % (cpt_colours[ftypes[name]]), width='0.2p')
        p.path(corners2gmt(corners[:, :2]), is_file=False)

    for x in p, q:
        x.finalise()
        x.png(dpi=300, background='white')
        call(['xdg-open', '%s.png' % (os.path.splitext(x.psf.name)[0])])
    raw_input('press return to continue...')
    rmtree(tmp)


def test_hypo_distribution(srf_dirs, out_dir):
    if os.path.isdir(out_dir):
        rmtree(out_dir)
    os.makedirs(out_dir)

    # normal distribution
    n_x = np.linspace(0, 1, 100)
    n_y = norm.cdf(n_x, 0.5, 0.25)
    # weibel distribution
    w_x = np.random.weibull(3.353, size=100000) * 0.612
    w_x.sort()
    w_y = np.arange(w_x.size) / (w_x.size - 1.0)

    for d in srf_dirs:
        n = os.path.basename(os.path.dirname(d))
        r = glob(os.path.join(d, '*.info'))
        shypo = np.zeros(len(r) + 2)
        dhypo = np.zeros(len(r) + 2)
        shypo[-1] = 1
        dhypo[-1] = 1
        for i, f in enumerate(r):
            with h5open(f, 'r') as h:
                shypo[i + 1] = h.attrs['shypo'][0] / sum(h.attrs['length'])
                dhypo[i + 1] = h.attrs['dhypo'][0] / h.attrs['width'][0]
        shypo.sort()
        dhypo.sort()
        y = np.arange(shypo.size) / (shypo.size - 1.0)

        fig = plt.figure(figsize = (8, 5), dpi = 150)
        plt.plot(n_x, n_y, label="Theoretical (along strike)")
        plt.plot(w_x, w_y, label="Theoretical (along dip)")
        # where = 'pre' or 'post'
        plt.step(shypo, y, label="Realizations (along strike)")
        plt.step(dhypo, y, label="Realizations (along dip)")
        plt.plot(dhypo[1:-1], y[1:-1], marker='.', markersize=2, linestyle='None')
        plt.plot(shypo[1:-1], y[1:-1], marker='.', markersize=2, linestyle='None')

        plt.legend(loc='best')
        plt.title('%s Hypocentre Location Distribution' % (n))
        plt.ylabel('CDF')
        plt.xlabel('along fault')
        plt.gca().set_ylim([0, 1])
        plt.gca().set_xlim([0, 1])
        plt.savefig(os.path.join(out_dir, n))
        plt.close()


parser = ArgumentParser()
parser.add_argument("nhm_srf_dir", help="NHM SRF directory to test")
parser.add_argument("--nhm-file", help="NHM fault list file", default=os.path.join(os.path.dirname(__file__), 'NZ_FLTmodel_2010.txt'))
parser.add_argument("--out-dir", help="Folder to place outputs in", default='./cstests')
parser.add_argument("-n", "--nproc", help="number of processes to use", type=int, default=1)
args = parser.parse_args()

assert os.path.isdir(args.nhm_srf_dir)
if not os.path.isdir(args.out_dir):
    os.makedirs(args.out_dir)

srf_dirs = glob(os.path.join(os.path.abspath(args.nhm_srf_dir), '*', 'Srf'))
#info_files = glob(os.path.join(os.path.abspath(args.nhm_srf_dir), '*', 'Srf', '*.info'))

#test_mw_vs_area(info_files)
#test_mw_vs_nrup(srf_dirs)
#test_seismogenic_depth(srf_dirs, args.nhm_file)
#test_spacial_srf(srf_dirs, args.nhm_file)
if args.nproc == 1:
    test_hypo_distribution(srf_dirs, os.path.join(args.out_dir, 'distributions'))
else:
    p = Pool(args.nproc)
    #p.map(test_hypo_distribution, srf_dirs)

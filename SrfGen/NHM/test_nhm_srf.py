#!/usr/bin/env python

from argparse import ArgumentParser
from glob import glob
import os

from h5py import File as h5open
import matplotlib.pyplot as plt
import numpy as np

from createSRF import leonard

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

parser = ArgumentParser()
parser.add_argument("nhm_srf_dir", help="NHM SRF directory to test")
args = parser.parse_args()

assert os.path.isdir(args.nhm_srf_dir)

srf_dirs = glob(os.path.join(os.path.abspath(args.nhm_srf_dir), '*', 'Srf'))
info_files = glob(os.path.join(os.path.abspath(args.nhm_srf_dir), '*', 'Srf', '*.info'))
test_mw_vs_area(info_files)
test_mw_vs_nrup(srf_dirs)

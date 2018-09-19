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
    fig = plt.figure(figsize = (16 / 3., 9 / 3.), dpi = 150)
    min_x = 10
    max_x = 0
    if len(ds_mw):
        min_x = min(min(ds_mw), min_x)
        max_x = max(max(ds_mw), max_x)
    if len(ss_mw):
        min_x = min(min(ss_mw), min_x)
        max_x = max(max(ss_mw), max_x)
    x = np.linspace(min_x, max_x, 100)

    if len(ds_mw):
        plt.plot(x, 10**(x-4), label='Leonard (dip slip)')
        plt.plot(ds_mw, ds_area, label='SRF (dip slip)', marker='x', linestyle='None')
    if len(ss_mw):
        plt.plot(x, 10**(x-3.99), label='Leonard (strike slip)')
        plt.plot(ss_mw, ss_area, label='SRF (strike slip)', marker='x', linestyle='None')
    plt.legend(loc='best')

    plt.gca().set_yscale('log')
    plt.show()
    plt.close()

parser = ArgumentParser()
parser.add_argument("nhm_srf_dir", help="NHM SRF directory to test")
args = parser.parse_args()

assert os.path.isdir(args.nhm_srf_dir)

info_files = glob(os.path.join(os.path.abspath(args.nhm_srf_dir), '*', 'Srf', '*.info'))
test_mw_vs_area(info_files)

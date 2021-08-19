#!/usr/bin/env python
import argparse
from glob import glob
import json
import math
import os
from typing import Dict, List

import pandas as pd
from os import path
from shutil import copy, rmtree
from tempfile import mkdtemp

from h5py import File as h5open
import matplotlib
from matplotlib.lines import Line2D
from qcore.nhm import load_nhm, NHMFault
from qcore.simulation_structure import get_fault_from_realisation

from srf_generation.pre_processing_common import load_realisation_file_as_dict
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    a_to_mw_leonard,
    mw_to_a_leonard,
    mw_sigma_leonard,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm, kstest
import numpy as np
import yaml

from qcore.geo import ll_shift, ll_dist
from qcore import gmt, srf, simulation_structure
from qcore.constants import VM_PARAMS_FILE_NAME


DS_RAKE = 90
SS_RAKE = 180


def fig_init():
    return plt.figure(figsize=(8, 5), dpi=150)


def format_and_save_plot(title=None, xlabel=None, ylabel=None, filepath=None):
    plt.legend(loc="best")
    if title is not None:
        plt.title(title)
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if filepath is not None:
        plt.savefig(filepath)
        plt.close()


def plot_parameter_histogram(
    param: str, realisations: Dict[str, pd.DataFrame], out_dir: str
):

    all_param = []
    mags = []

    for fault, table in realisations.items():
        all_param.extend(table[param])
        if param is not "magnitude":
            mags.extend(table["magnitude"])

    fig_init()
    plt.hist(all_param)
    format_and_save_plot(
        title=f"{param} distribution. Mean: {np.mean(all_param)}",
        xlabel=f"{param}",
        ylabel=f"n",
        filepath=path.join(out_dir, f"{param}_hist"),
    )
    if param is not "magnitude":
        fig_init()
        plt.plot(all_param, mags)
        format_and_save_plot(
            title=f"{param} vs magnitude",
            xlabel=f"{param}",
            ylabel="Magnitude",
            filepath=path.join(out_dir, f"{param}_vs_mag"),
        )


###
### SRF PLOTS
###


def plot_area_mag(
    faults: Dict[str, pd.DataFrame],
    median_faults: Dict[str, pd.DataFrame],
    out_dir: str,
):
    # load mw and areas from srf infos
    ds_mw = []
    ds_mw_mean = []
    ds_area = []
    ds_area_mean = []
    ss_mw = []
    ss_mw_mean = []
    ss_area = []
    ss_area_mean = []
    for fault, df in faults.items():
        area = df["length"] * df["width_subfault_0"]
        # if dip slip else strike slip
        if round(median_faults[fault]["rake"][0] % 360 / 90.0) % 2:
            ds_mw.extend(df["magnitude"])
            ds_area.extend(area)
            ds_mw_mean.extend(median_faults[fault]["magnitude"])
            ds_area_mean.extend(
                median_faults[fault]["length"]
                * median_faults[fault]["width_subfault_0"]
            )
        else:
            ss_mw.extend(df["magnitude"])
            ss_area.extend(area)
            ss_mw_mean.extend(median_faults[fault]["magnitude"])
            ss_area_mean.extend(
                median_faults[fault]["length"]
                * median_faults[fault]["width_subfault_0"]
            )

    # plot
    fig_init()
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
        plt.semilogx(mw_to_a_leonard(y, 90), y, label="Leonard (dip slip) - 1 s.d.")
        plt.semilogx(mw_to_a_leonard(y, 90), y, label="Leonard (dip slip)")
        plt.semilogx(
            mw_to_a_leonard(y, 90) + mw_sigma_leonard(90),
            y,
            label="Leonard (dip slip) + 1 s.d.",
        )
        plt.semilogx(
            ds_area, ds_mw, label="SRF (dip slip)", marker="x", linestyle="None"
        )
        plt.semilogx(
            ds_area_mean,
            ds_mw_mean,
            label="SRF (dip slip) mean",
            marker="x",
            linestyle="None",
        )
    if len(ss_mw):
        plt.semilogx(mw_to_a_leonard(y, 180), y, label="Leonard (strike slip)")
        plt.semilogx(
            ss_area, ss_mw, label="SRF (strike slip)", marker="x", linestyle="None"
        )
        plt.semilogx(
            ss_area_mean,
            ss_mw_mean,
            label="SRF (strike slip) mean",
            marker="x",
            linestyle="None",
        )

    format_and_save_plot(
        title="Area vs Magnitude",
        xlabel="Area (sq. km)",
        ylabel="Magnitude from SRF info",
        filepath=os.path.join(out_dir, "area_vs_mag"),
    )

    # part 2 : mw cs vs mw leonard
    fig_init()
    plt.plot(y, y - mw_sigma_leonard(90), label="Target - 1 s.d.")
    plt.plot(y, y, label="Target")
    plt.plot(y, y + mw_sigma_leonard(90), label="Target + 1 s.d.")
    if len(ds_mw):
        plt.plot(
            ds_mw,
            a_to_mw_leonard(ds_area, 4, 3.99, 90),
            label="Actual (dip slip)",
            marker="x",
            linestyle="None",
        )
        plt.plot(
            ds_mw_mean,
            a_to_mw_leonard(ds_area_mean, 4, 3.99, 90),
            label="Actual (dip slip) (unperturbated)",
            marker="x",
            linestyle="None",
        )
    if len(ss_mw):
        plt.plot(
            ss_mw,
            a_to_mw_leonard(ss_area, 4, 3.99, 180),
            label="Actual (strike slip)",
            marker="x",
            linestyle="None",
        )
        plt.plot(
            ss_mw_mean,
            a_to_mw_leonard(ss_area_mean, 4, 3.99, 180),
            label="Actual (strike slip) (unperturbated)",
            marker="x",
            linestyle="None",
        )
    format_and_save_plot(
        title="Mw SRF vs Mw Leonard",
        xlabel="Mw from SRF info",
        ylabel="Mw based on area, based on Leonard",
        filepath=os.path.join(out_dir, "mag_vs_leonard"),
    )


def plot_area_nhm(
    median_faults: Dict[str, pd.DataFrame], nhm_file: Dict[str, NHMFault], out_dir: str
):
    srf_area = []
    fault_name = []
    for fault, df in median_faults.items():
        area = df["length"] * df["width_subfault_0"]
        srf_area.append(area[0])
        fault_name.append(fault)

    nhm_area = []
    for fault in fault_name:
        nhm_info: NHMFault = nhm_file[fault]
        dlen = sum(
            [
                ll_dist(
                    nhm_info.trace[i][0],
                    nhm_info.trace[i][1],
                    nhm_info.trace[i + 1][0],
                    nhm_info.trace[i + 1][1],
                )
                for i in range(len(nhm_info.trace) - 1)
            ]
        )
        dwid = (nhm_info.dbottom - nhm_info.dtop) / math.sin(math.radians(nhm_info.dip))
        nhm_area.append(dlen * dwid)

    # plot
    fig_init()
    plt.plot(
        [min(nhm_area), max(nhm_area)], [min(nhm_area), max(nhm_area)], label="Target"
    )
    plt.plot(nhm_area, srf_area, label="Actual", marker="x", linestyle="None")
    format_and_save_plot(
        title="Area NHM vs Area SRF",
        xlabel="Area from NHM (km sq.)",
        ylabel="Area from SRF (km sq.)",
        filepath=os.path.join(out_dir, "area_vs_nhm"),
    )

    # text
    with open(os.path.join(out_dir, "area_vs_nhm.txt"), "w") as t:
        t.write("name nhm srf\n")
        for i, name in enumerate(fault_name):
            t.write(f"{name} {nhm_area[i]:4.2} {srf_area[i]:4.2}\n")


def plot_mag_nhm(
    faults: Dict[str, pd.DataFrame],
    median_faults: Dict[str, pd.DataFrame],
    nhm_file: Dict[str, NHMFault],
    out_dir: str,
):
    # part 3 : mw cs vs mw nhm

    srf_mw = np.zeros((len(faults), 50)) * np.nan
    srf_mw_median = []
    nhm_mw = []

    for i, (fault, df) in enumerate(faults.items()):
        srf_mw[i][: faults[fault].shape[0]] = df["magnitude"]
        srf_mw_median.extend(median_faults[fault]["magnitude"])
        nhm_mw.append(nhm_file[fault].mw)

    # plot
    fig_init()
    plt.plot(
        [
            np.nanmin(
                [np.nanmin(flattened) for flattened in srf_mw] + srf_mw_median + nhm_mw
            ),
            np.nanmax(
                [np.nanmax(flattened) for flattened in srf_mw] + srf_mw_median + nhm_mw
            ),
        ],
        [
            np.nanmin(
                [np.nanmin(flattened) for flattened in srf_mw] + srf_mw_median + nhm_mw
            ),
            np.nanmax(
                [np.nanmax(flattened) for flattened in srf_mw] + srf_mw_median + nhm_mw
            ),
        ],
        label="Target",
    )
    points: List[Line2D] = plt.plot(
        nhm_mw, srf_mw, "bx", label="Actual", linestyle="None"
    )
    for p in points[1:]:
        p.set_label(None)
    plt.plot(nhm_mw, srf_mw_median, label="Actual median", marker="x", linestyle="None")

    format_and_save_plot(
        title="Mw NHM vs Mw SRF info3",
        xlabel="Mw from NHM",
        ylabel="Mw from SRF info",
        filepath=os.path.join(out_dir, "mag_vs_nhm"),
    )


def plot_mag_nrup(
    faults: Dict[str, pd.DataFrame],
    median_faults: Dict[str, pd.DataFrame],
    out_dir: str,
):
    rel_count = []
    srf_mw_median = []

    for fault, df in faults.items():
        rel_count.append(df.shape[0])
        srf_mw_median.extend(median_faults[fault]["magnitude"])

    # plot
    x = [min(srf_mw_median), max(srf_mw_median)]
    if x[0] < 6 < x[-1]:
        x.insert(1, 6)
    if x[0] < 8 < x[-1]:
        x.insert(-1, 8)

    fig_init()
    plt.plot(
        x, np.maximum(np.minimum((np.array(x) * 20 - 110), 50), 10), label="Target"
    )
    plt.plot(srf_mw_median, rel_count, label="Actual", marker="x", linestyle="None")

    format_and_save_plot(
        title="Magnitude vs Number of Realizations",
        xlabel="Magnitude",
        ylabel="Number of Realizations",
        filepath=os.path.join(out_dir, "mag_vs_nrup"),
    )


def plot_dbottom(info_files: List[str], nhm_data: Dict[str, NHMFault], out_dir: str):

    srf_depths = []
    fault_names = []
    for info_file in info_files:
        with h5open(info_file, "r") as i:
            srf_depths.append(np.average(i.attrs["dbottom"]))
        fault_names.append(os.path.basename(info_file).split(".")[0])

    # grab depths from nhm
    nhm_depths = np.zeros(len(fault_names))
    for fault in fault_names:
        nhm_info = nhm_data[get_fault_from_realisation(fault)]
        nhm_depths[fault_names.index(fault)] = nhm_info.dbottom

    # plot
    fig_init()
    max_x = max(max(srf_depths), max(nhm_depths))
    plt.plot([0, max_x], [0, max_x], label="NHM")
    plt.plot([0, max_x], [3, max_x + 3], label="NHM + 3")
    plt.plot(
        nhm_depths,
        np.array(srf_depths),
        label="SRF Depths",
        marker="x",
        linestyle="None",
    )

    format_and_save_plot(
        title="Seismogenic Depth from NHM and SRF",
        xlabel="NHM Depth (km)",
        ylabel="SRF Depth (km)",
        filepath=os.path.join(out_dir, "dbottom"),
    )


def nhm2corners(nhm_data: Dict[str, NHMFault], names: List[str]):
    corners_sub = []
    corners_shal = []
    ftypes = {}
    ex_names = []
    for n in nhm_data.keys():
        if n in names:
            ftypes[n] = 2 + (not nhm_data[n].fault_type.startswith("SUBDUCTION"))
        else:
            fault_data = nhm_data[n]
            ex_names.append(n)
            pwid = (fault_data.dbottom - fault_data.dtop) / np.tan(
                np.radians(fault_data.dip)
            )

            ftype = int(not nhm_data[n].fault_type.startswith("SUBDUCTION"))
            corners = np.zeros((len(fault_data.trace) - 1, 4, 2))
            trace = fault_data.trace
            dip_dir = fault_data.dip_dir
            for i in range(len(corners)):
                corners[i, :2] = trace[i : i + 2]
                corners[i, 2] = ll_shift(
                    corners[i, 1, 1], corners[i, 1, 0], pwid, dip_dir
                )[::-1]
                corners[i, 3] = ll_shift(
                    corners[i, 0, 1], corners[i, 0, 0], pwid, dip_dir
                )[::-1]
            if ftype % 2:
                corners_shal.append(corners)
            else:
                corners_sub.append(corners)

    return ftypes, corners_shal, corners_sub, ex_names


def plot_srf_trace_dbottom(info_files, out_dir):
    def trace2gmt(trace, depth):
        return "> -Z%s\n" % (depth) + "\n".join(
            [" ".join(map(str, ll)) for ll in trace]
        )

    tmp = mkdtemp()
    cpt = os.path.join(tmp, "cpt.cpt")
    depth_min = 5
    depth_max = 30
    gmt.makecpt(
        "rainbow",
        cpt,
        depth_min,
        depth_max,
        inc=(depth_max - depth_min) / 5.0,
        continuous=False,
    )

    p = gmt.GMTPlot(os.path.join(tmp, "trace.ps"))
    p.spacial("M", region=gmt.nz_region, sizing=10, x_shift=1, y_shift=1)
    p.basemap(topo=None, road=None, highway=None, land="white", water="white", res="f")
    p.ticks(major="2d")
    p.text(
        sum(gmt.nz_region[:2]) / 2.0,
        gmt.nz_region[3],
        "Fault Depths",
        size="26p",
        dy=0.3,
    )

    # fault traces, colour based on bottom depth
    for info_file in info_files:
        with h5open(info_file, "r") as h:
            corners = h.attrs["corners"]
            dbottom = h.attrs["dbottom"]
            p.path(
                "\n".join(
                    [
                        trace2gmt(corners[i, :2], depth=dbottom[i])
                        for i in range(len(corners))
                    ]
                ),
                is_file=False,
                colour="40/40/40",
                cpt=cpt,
                width="1.2p",
            )

    p.cpt_scale(
        "C",
        "R",
        cpt,
        label="depth to bottom (km)",
        pos="rel_out",
        horiz=False,
        arrow_f=False,
        length=5,
        thickness=0.3,
        dx=0.3,
        categorical=True,
        intervals=True,
        gap="0.3",
        zmax=depth_max,
    )

    p.finalise()
    p.png(dpi=300, background="white", out_dir=out_dir)


def plot_srf_nhm(info_files, nhm_file, out_dir):
    tmp = mkdtemp()
    cpt = os.path.join(tmp, "types.cpt")
    cpt_labels = {
        "0": ";subduction excluded\n",
        "1": ";shallow excluded\n",
        "2": ";subduction included\n",
        "3": ";shallow included\n",
    }
    gmt.makecpt("rainbow", cpt, 0, 4, inc=1, continuous=False)
    cpt_colours = []
    with open(cpt, "r") as c:
        cpt_lines = c.readlines()
    for i in range(len(cpt_lines)):
        start = cpt_lines[i].split()[0]
        if start in cpt_labels:
            cpt_colours.append(cpt_lines[i].split()[1])
            cpt_lines[i] = "%s%s" % (cpt_lines[i].strip(), cpt_labels[start])
    with open(cpt, "w") as c:
        c.write("".join(cpt_lines))

    p = gmt.GMTPlot(os.path.join(tmp, "srfs.ps"))
    p.spacial("M", region=gmt.nz_region, sizing=10, x_shift=1, y_shift=1)
    p.basemap(
        topo=None, road=None, highway=None, land="lightgray", water="white", res="f"
    )
    p.ticks(major="2d")
    p.leave()
    copy(p.psf.name, os.path.join(tmp, "srfs_ex.ps"))
    q = gmt.GMTPlot(os.path.join(tmp, "srfs_ex.ps"), append=True, reset=False)
    p.enter()

    title = "SRF and rest of NHM faults"
    zmax = "NaN"
    length = 8 * 0.3
    for x in p, q:
        x.text(
            sum(gmt.nz_region[:2]) / 2.0, gmt.nz_region[3], title, size="26p", dy=0.3
        )
        x.cpt_scale(
            "C",
            "R",
            cpt,
            pos="rel_out",
            horiz=False,
            arrow_f=False,
            length=length,
            thickness=0.3,
            dx=0.3,
            categorical=True,
            gap="0.3",
            zmax=zmax,
        )
        title = "Rest of NHM faults"
        zmax = 2
        length = 4 * 0.3

    def corners2gmt(corners):
        return "\n>\n".join(
            ["\n".join([" ".join(map(str, ll)) for ll in plane]) for plane in corners]
        )

    # nhm faults
    faults = sorted(
        list(
            {
                get_fault_from_realisation(path.basename(f_name).split(".")[0])
                for f_name in info_files
            }
        )
    )
    ftypes, nhm_cnrs_shal, nhm_cnrs_sub, ex_names = nhm2corners(nhm_file, faults)
    for x in p, q:
        x.path(
            "\n>\n".join(map(corners2gmt, nhm_cnrs_shal)),
            is_file=False,
            colour="40/40/40",
            split="-",
            fill="%s@80" % (cpt_colours[1]),
            close=True,
            width="0.2p",
        )
        x.path(
            "\n>\n".join(map(corners2gmt, nhm_cnrs_sub)),
            is_file=False,
            colour="40/40/40",
            split="-",
            fill="%s@80" % (cpt_colours[0]),
            close=True,
            width="0.2p",
        )
        x.path(
            "\n>\n".join([corners2gmt(c[:, :2]) for c in nhm_cnrs_shal]),
            is_file=False,
            colour="40/40/40",
        )
        x.path(
            "\n>\n".join([corners2gmt(c[:, :2]) for c in nhm_cnrs_sub]),
            is_file=False,
            colour="40/40/40",
        )

    # srf files
    for i, f in enumerate(info_files):
        name = get_fault_from_realisation(path.basename(f).split(".")[0])
        with h5open(f, "r") as h:
            corners = h.attrs["corners"]
        p.path(
            corners2gmt(corners),
            is_file=False,
            close=True,
            split="-",
            fill="%s@80" % (cpt_colours[ftypes[name]]),
            width="0.2p",
        )
        p.path(corners2gmt(corners[:, :2]), is_file=False)

    for x in p, q:
        x.finalise()
        x.png(dpi=300, background="white", out_dir=out_dir)
    rmtree(tmp)

    # save excluded fault names
    with open(os.path.join(out_dir, "excluded_faults.txt"), "w") as e:
        e.write("\n".join(ex_names))
        e.write("\n")


def plot_hypo_dist(
    faults: Dict[str, pd.DataFrame],
    median_faults: Dict[str, pd.DataFrame],
    out_dir: str,
):
    if os.path.isdir(out_dir):
        rmtree(out_dir)
    os.makedirs(out_dir)

    # normal distribution
    n_x = np.linspace(0, 1, 100)
    n_y = norm.cdf(n_x, 0.5, 0.25)
    n_call = lambda x: norm.cdf(x, 0.5, 0.25)
    # weibull distribution
    w_x = np.random.weibull(3.353, size=100_000) * 0.612
    w_x.sort()
    w_x = np.asarray(w_x)[(0 <= np.asarray(w_x)) * (np.asarray(w_x) <= 1)]
    w_y = np.arange(w_x.size) / (w_x.size - 1.0)
    w_call = lambda x: w_y[[np.argmin(np.abs(v - w_x)) for v in x]]

    # kstest
    names = list(faults.keys())
    mw = []
    p_s = []
    p_d = []
    for fault_name in names:
        fault_data = faults[fault_name]

        shypo = sorted(
            (
                fault_data["shypo"]
                / sum(
                    [
                        fault_data[f"length_subfault_{i}"]
                        for i in range(fault_data["plane_count"][0])
                    ]
                )
            ).to_numpy("float")
            + 0.5
        )
        p_s.append(kstest(shypo, n_call).pvalue)
        dhypo = sorted(
            (fault_data["dhypo"] / fault_data["width_subfault_0"]).to_numpy("float")
        )
        p_d.append(kstest(dhypo, w_call).pvalue)
        mw.extend(median_faults[fault_name]["magnitude"])

        y = np.arange(len(shypo)) / (len(shypo) - 1.0)

        fig_init()
        plt.plot(
            n_x,
            n_y,
            label="Theoretical (along strike)",
            color="blue",
            linestyle="dashed",
        )
        plt.step(shypo, y, label="Realizations (along strike)", color="blue")
        plt.plot(
            w_x, w_y, label="Theoretical (along dip)", color="red", linestyle="dashed"
        )
        plt.step(dhypo, y, label="Realizations (along dip)", color="red")
        plt.plot(
            dhypo[1:-1],
            y[1:-1],
            marker=".",
            markersize=2,
            linestyle="None",
            color="orange",
        )
        plt.plot(
            shypo[1:-1],
            y[1:-1],
            marker=".",
            markersize=2,
            linestyle="None",
            color="green",
        )

        plt.gca().set_ylim([0, 1])
        plt.gca().set_xlim([0, 1])

        format_and_save_plot(
            title=f"{fault_name} Hypocentre Location Distribution ({faults[fault_name].shape[0]} realizations)",
            xlabel="along fault",
            ylabel="CDF",
            filepath=os.path.join(out_dir, fault_name),
        )

    # overall plot
    fig_init()
    plt.plot([min(mw), max(mw)], [0.05, 0.05])
    plt.plot(
        mw,
        p_s,
        label="Num p (along strike) < 0.05 = %d" % (sum(np.array(p_s) < 0.05)),
        marker="x",
        linestyle="None",
    )
    plt.plot(
        mw,
        p_d,
        label="Num p (along dip) < 0.05 = %d" % (sum(np.array(p_d) < 0.05)),
        marker="x",
        linestyle="None",
    )

    plt.gca().set_yscale("log")
    # plt.gca().set_ylim([0, 1])
    format_and_save_plot(
        title=f"Hypocentre Location Distribution ({len(mw)} simulations)",
        xlabel="magnitude",
        ylabel="p",
        filepath=os.path.join(out_dir, "..", "hypo_dist"),
    )

    # overall histogram strike
    fig_init()
    plt.hist(p_s, bins=1 + int(round(max(p_s) / 0.05)), label="along strike")
    format_and_save_plot(
        title=f"Hypocentre along Strike Distribution Probability of Fit ({len(mw)} simulations)",
        xlabel="probability of fitting normal distribution",
        ylabel="n",
        filepath=os.path.join(out_dir, "..", "hypo_dist_strike"),
    )

    # overall histogram dip
    fig_init()
    plt.hist(p_d, bins=1 + int(round(max(p_d) / 0.05)), label="along dip")
    format_and_save_plot(
        title=f"Hypocentre along Dip Distribution Probability of Fit ({len(mw)} simulations)",
        xlabel="probability of fitting weibull distribution",
        ylabel="n",
        filepath=os.path.join(out_dir, "..", "hypo_dist_dip"),
    )

    # as text
    with open(os.path.join(out_dir, "..", "hypo_dist.txt"), "w") as d:
        d.write("name strike-p dip-p mw\n")
        for i in range(len(names)):
            d.write("%s %s %s %s\n" % (names[i], p_s[i], p_d[i], mw[i]))


def plot_srf_error(info_files: List[str], out_dir: str):
    # strike and dip error assuming 0.1km subfault spacing
    error_s = []
    error_d = []
    for i in info_files:
        with h5open(i, "r") as h:
            error_s.append(
                (sum(h.attrs["nstrike"]) * 0.1 - sum(h.attrs["length"]))
                / sum(h.attrs["length"])
            )
            error_d.append(
                (h.attrs["ndip"][0] * 0.1 - h.attrs["width"][0]) / h.attrs["width"][0]
            )

    fig_init()
    plt.plot(np.arange(len(error_s)), error_s, label="", marker="x", linestyle="None")
    format_and_save_plot(
        title="Difference in Length / Length (nstrike * 0.1km)",
        xlabel="sources",
        ylabel="SRF length / length difference",
        filepath=os.path.join(out_dir, "srf_hh_length"),
    )
    #
    fig_init()
    plt.plot(np.arange(len(error_d)), error_d, label="", marker="x", linestyle="None")
    format_and_save_plot(
        title="Difference in Width / Width (ndip * 0.1km)",
        xlabel="sources",
        ylabel="SRF width / width difference",
        filepath=os.path.join(out_dir, "srf_hh_width"),
    )


###
### VM PLOTS
###


def plot_vm(vm_dirs, out_dir):
    # load vm corners
    vm_corners = [os.path.join(d, "VeloModCorners.txt") for d in vm_dirs]
    vm_corners = "\n>\n".join(
        [
            "\n".join(
                [
                    " ".join(map(str, v))
                    for v in np.loadtxt(c, skiprows=2, dtype=np.float32).tolist()
                ]
            )
            for c in vm_corners
        ]
    )

    tmp = mkdtemp()
    p = gmt.GMTPlot(os.path.join(tmp, "corners.ps"))
    p.spacial("M", region=gmt.nz_region, sizing=10, x_shift=1, y_shift=1)
    p.basemap(
        topo=None, road=None, highway=None, land="lightgray", water="white", res="f"
    )
    p.ticks(major="2d")
    p.text(
        sum(gmt.nz_region[:2]) / 2.0, gmt.nz_region[3], "VM Corners", size="26p", dy=0.3
    )

    # plot velocity model corners
    p.path(vm_corners, is_file=False, close=True, width="0.5p", split="-")

    p.finalise()
    p.png(dpi=300, background="white", out_dir=out_dir)
    rmtree(tmp)


def plot_duration_mag(vm_dirs, out_dir):
    durations = []
    magnitudes = []
    for j in vm_dirs:
        with open(path.join(j, "vm_params.yaml"), "r") as jo:
            d = yaml.load(jo)
        magnitudes.append(d["mag"])
        durations.append(d["sim_duration"])

    fig_init()
    plt.plot(
        magnitudes, durations, label="Velocity Models", marker="x", linestyle="None"
    )

    format_and_save_plot(
        title="Magnitude vs Simulation Duration",
        xlabel="Magnitude",
        ylabel="Simulation Duration (s)",
        filepath=os.path.join(out_dir, "mag_vs_duration"),
    )


def plot_srf_subfault_size(srf_files, out_dir):
    srf_subfault_size = []
    for i, srf_file in enumerate(srf_files):
        i_dx, i_dy = srf.srf_dxy(srf_file)
        srf_subfault_size.append((i, i_dx, i_dy))

    x, dx, dy = zip(*srf_subfault_size)

    fig_init()
    plt.plot(x, dx, marker="x", linestyle="None")
    format_and_save_plot(
        title="SRF Subfault Size (dx)",
        xlabel="sources",
        ylabel="dx (km)",
        filepath=os.path.join(out_dir, "srf_dx_subfault_size"),
    )

    fig_init()
    plt.plot(x, dy, marker="x", linestyle="None")
    format_and_save_plot(
        title="SRF Subfault Size (dy)",
        xlabel="sources",
        ylabel="dy (km)",
        filepath=os.path.join(out_dir, "srf_dy_subfault_size"),
    )


def plot_stoch_subfault_size(stoch_files, out_dir):
    stoch_subfault_size = []
    for i, stoch_file in enumerate(stoch_files):
        stoch_header = srf.read_stoch_header(stoch_file)
        for plane in stoch_header["plane_headers"]:
            stoch_subfault_size.append((i, plane["dx"], plane["dy"]))
    x, dx, dy = zip(*stoch_subfault_size)

    fig_init()
    plt.plot(x, dx, marker="x", linestyle="None")
    format_and_save_plot(
        title="Stoch Subfault Size (dx)",
        xlabel="sources",
        ylabel="dx (km)",
        filepath=os.path.join(out_dir, "stoch_dx_subfault_size"),
    )

    fig_init()
    plt.plot(x, dy, marker="x", linestyle="None")
    format_and_save_plot(
        title="Stoch Subfault Size (dy)",
        xlabel="sources",
        ylabel="dy (km)",
        filepath=os.path.join(out_dir, "stoch_dy_subfault_size"),
    )


def plot_size_mag(vm_dirs, out_dir):
    sizes = []
    magnitudes = []
    for d in vm_dirs:
        sizes.append(os.stat(os.path.join(d, "vs3dfile.s")).st_size / 1_000_000.0)
        with open(os.path.join(d, VM_PARAMS_FILE_NAME), "r") as jo:
            magnitudes.append(yaml.load(jo)["mag"])

    fig_init()
    plt.plot(
        magnitudes,
        sizes,
        label="Velocity Model Components",
        marker="x",
        linestyle="None",
    )

    format_and_save_plot(
        title="Magnitude vs VM File Sizes",
        xlabel="Magnitude",
        ylabel="File Sizes (megabytes)",
        filepath=os.path.join(out_dir, "mag_vs_size"),
    )


def plot_vm_params(vm_dirs, out_dir):
    params = {"hh": [], "flo": []}
    for j in vm_dirs:
        with open(path.join(j, "vm_params.yaml"), "r") as jo:
            d = json.load(jo)
            params["hh"].append(d["hh"])
            params["flo"].append(d["flo"])
    x = np.arange(len(vm_dirs))

    for parameter in params:
        fig_init()
        plt.plot(x, params[parameter], label=parameter, marker="x", linestyle="None")
        plt.gca().get_xaxis().set_visible(False)
        format_and_save_plot(
            title=f"VM Parameter ({parameter})",
            ylabel=parameter,
            filepath=os.path.join(out_dir, f"vm_{parameter}"),
        )


def load_rel_files(rel_files):
    faults = {}
    current_fault_name = ""
    current_fault_rels = {}

    for file_name in sorted(rel_files):
        realisation = load_realisation_file_as_dict(file_name)
        name_parts = realisation["name"].split("_")
        if name_parts[0] != current_fault_name:
            if len(current_fault_rels) > 0:
                faults[str(current_fault_name)] = pd.DataFrame(
                    current_fault_rels
                ).transpose()
                current_fault_rels = {}
            current_fault_name = name_parts[0]
        current_fault_rels[name_parts[1] if len(name_parts) > 1 else 0] = realisation
    faults[str(current_fault_name)] = pd.DataFrame(current_fault_rels).transpose()

    return faults


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("nhm_file", type=path.abspath)
    parser.add_argument("--params_to_plot", nargs="+", type=str, default=[])
    parser.add_argument(
        "--cybershake_root", type=path.abspath, default=path.abspath(".")
    )
    parser.add_argument(
        "--out-dir",
        help="Folder to place outputs in",
        default=path.abspath("./rel_validation_plots"),
    )
    return parser.parse_args()


def main():
    args = load_args()

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    nhm_data = load_nhm(args.nhm_file)

    stoch_files = sorted(
        glob(simulation_structure.get_stoch_path(args.cybershake_root, "*"))
    )

    srf_files = sorted(
        glob(simulation_structure.get_srf_path(args.cybershake_root, "*"))
    )
    info_files = list(map(lambda x: x.replace(".srf", ".info"), srf_files))
    rel_files = list(map(lambda x: x.replace(".srf", ".csv"), srf_files))
    faults = load_rel_files(rel_files)

    median_rel_files = glob(
        path.join(
            simulation_structure.get_sources_dir(args.cybershake_root), "*", "*.csv"
        )
    )
    median_faults = load_rel_files(median_rel_files)

    for param in args.params_to_plot:
        plot_parameter_histogram(param, faults, args.out_dir)

    plot_area_mag(faults, median_faults, args.out_dir)
    plot_area_nhm(median_faults, nhm_data, args.out_dir)
    plot_mag_nhm(faults, median_faults, nhm_data, args.out_dir)
    plot_mag_nrup(faults, median_faults, args.out_dir)

    plot_dbottom(info_files, nhm_data, args.out_dir)
    plot_srf_trace_dbottom(info_files, args.out_dir)
    plot_srf_nhm(info_files, nhm_data, args.out_dir)

    plot_hypo_dist(faults, median_faults, os.path.join(args.out_dir, "hypo_dist"))
    plot_srf_error(info_files, args.out_dir)
    plot_stoch_subfault_size(stoch_files, args.out_dir)
    plot_srf_subfault_size(srf_files, args.out_dir)

    vm_dirs = sorted(
        glob(simulation_structure.get_fault_VM_dir(args.cybershake_root, "*"))
    )

    # # run VM plots and checks
    if len(vm_dirs):

        plot_vm(vm_dirs, args.out_dir)
        plot_duration_mag(vm_dirs, args.out_dir)
        plot_size_mag(vm_dirs, args.out_dir)
        plot_vm_params(vm_dirs, args.out_dir)


if __name__ == "__main__":
    main()

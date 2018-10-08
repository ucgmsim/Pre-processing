#!/usr/bin/env python

from argparse import ArgumentParser
from glob import glob
import json
import math
import os
from shutil import copy, rmtree
from subprocess import call
from tempfile import mkdtemp

from h5py import File as h5open
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm, kstest
import numpy as np

from createSRF import leonard
from qcore.geo import ll_shift, ll_bearing, ll_dist
from qcore import gmt

def fig_init():
    return plt.figure(figsize=(8, 5), dpi=150)

###
### SRF PLOTS
###

def plot_area_mag(srf_infos, out_dir):
    # load mw and areas from srf infos
    ds_mw = []
    ds_area = []
    ss_mw = []
    ss_area = []
    for i in srf_infos:
        with h5open(i, "r") as h:
            area = sum(h.attrs["length"] * h.attrs["width"])
            # if dip slip else strike slip
            if round(h.attrs["rake"][0][0] % 360 / 90.0) % 2:
                ds_mw.append(h.attrs["mag"])
                ds_area.append(area)
            else:
                ss_mw.append(h.attrs["mag"])
                ss_area.append(area)

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
        plt.plot(10 ** (y - 4), y, label="Leonard (dip slip)")
        plt.plot(ds_area, ds_mw, label="SRF (dip slip)", marker="x", linestyle="None")
    if len(ss_mw):
        plt.plot(10 ** (y - 3.99), y, label="Leonard (strike slip)")
        plt.plot(
            ss_area, ss_mw, label="SRF (strike slip)", marker="x", linestyle="None"
        )

    plt.legend(loc="best")
    plt.gca().set_xscale("log")
    plt.title("Area vs Magnitude")
    plt.ylabel("Magnitude from SRF info")
    plt.xlabel("Area (sq. km)")
    plt.savefig(os.path.join(out_dir, "area_vs_mag"))
    plt.close()

    # part 2 : mw cs vs mw leonard
    fig_init()
    plt.plot(y, y, label="Target")
    if len(ds_mw):
        plt.plot(
            ds_mw,
            4 + np.log10(ds_area),
            label="Actual (dip slip)",
            marker="x",
            linestyle="None",
        )
    if len(ss_mw):
        plt.plot(
            ss_mw,
            3.99 + np.log10(ss_area),
            label="Actual (strike slip)",
            marker="x",
            linestyle="None",
        )
    plt.legend(loc="best")
    plt.title("Mw SRF vs Mw Leonard")
    plt.xlabel("Mw from SRF info")
    plt.ylabel("Mw based on area, based on Leonard")
    plt.savefig(os.path.join(out_dir, "mag_vs_leonard"))
    plt.close()


def plot_area_nhm(srf_dirs, nhm_file, out_dir):
    area_srf_info = []
    fault_names = []
    for d in srf_dirs:
        info = glob(os.path.join(d, "*.info"))[0]
        with h5open(info, "r") as i:
            area_srf_info.append(sum(i.attrs["length"]) * i.attrs["width"][0])
        fault_names.append(os.path.basename(os.path.dirname(d)))

    # grab area from nhm
    area_nhm = np.zeros(len(fault_names))
    with open(nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
        n_i = 0
        while min(area_nhm) <= 0:
            try:
                name = nhm[n_i].strip()
            except IndexError:
                print("Could not find all SRFs in NHM.")
                return
            n_pt = int(nhm[n_i + 11])
            if name in fault_names:
                strace = [
                    map(float, ll)
                    for ll in map(
                        str.split, map(str.strip, nhm[n_i + 12 : n_i + 12 + n_pt])
                    )
                ]
                slen = sum(
                    [
                        ll_dist(
                            strace[i][0],
                            strace[i][1],
                            strace[i + 1][0],
                            strace[i + 1][1],
                        )
                        for i in range(n_pt - 1)
                    ]
                )
                dbottom = float(nhm[n_i + 6].split()[0])
                if dbottom <= 12:
                    dbottom += 3
                dwid = (dbottom - float(nhm[n_i + 7].split()[0])) / math.sin(
                    math.radians(float(nhm[n_i + 3].split()[0]))
                )
                area_nhm[fault_names.index(name)] = slen * dwid
            n_i += 13 + n_pt

    # plot
    fig_init()
    plt.plot(
        [min(area_nhm), max(area_nhm)], [min(area_nhm), max(area_nhm)], label="Target"
    )
    plt.plot(area_nhm, area_srf_info, label="Actual", marker="x", linestyle="None")
    plt.legend(loc="best")
    plt.title("Area NHM vs Area SRF")
    plt.xlabel("Area from NHM (km sq.)")
    plt.ylabel("Area from SRF (km sq.)")
    plt.savefig(os.path.join(out_dir, "area_vs_nhm"))
    plt.close()

    # text
    with open(os.path.join(out_dir, "area_vs_nhm.txt"), "w") as t:
        t.write("name nhm srf\n")
        for i in range(len(fault_names)):
            t.write("%s %s %s\n" % (fault_names[i], area_nhm[i], area_srf_info[i]))


def plot_mag_nhm(srf_dirs, nhm_file, out_dir):
    # part 3 : mw cs vs mw nhm
    mw_srf_info = []
    fault_names = []
    for d in srf_dirs:
        info = glob(os.path.join(d, "*.info"))[0]
        with h5open(info, "r") as i:
            mw_srf_info.append(i.attrs["mag"])
        fault_names.append(os.path.basename(os.path.dirname(d)))

    # grab mw from nhm
    mw_nhm = np.zeros(len(fault_names))
    with open(nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
        n_i = 0
        while min(mw_nhm) <= 0:
            name = nhm[n_i].strip()
            if name in fault_names:
                mw_nhm[fault_names.index(name)] = float(nhm[n_i + 10].split()[0])
            elif name == "":
                print("Could not find all SRFs in NHM.")
                return
            n_i += 13 + int(nhm[n_i + 11])

    # plot
    fig_init()
    plt.plot([min(mw_nhm), max(mw_nhm)], [min(mw_nhm), max(mw_nhm)], label="Target")
    plt.plot(mw_nhm, mw_srf_info, label="Actual", marker="x", linestyle="None")
    plt.legend(loc="best")
    plt.title("Mw NHM vs Mw SRF info")
    plt.xlabel("Mw from NHM")
    plt.ylabel("Mw from SRF info")
    plt.savefig(os.path.join(out_dir, "mag_vs_nhm"))
    plt.close()


def plot_mag_nrup(srf_dirs, out_dir):
    nrup = []
    mw = []
    for d in srf_dirs:
        infos = glob(os.path.join(d, "*.info"))
        nrup.append(len(infos))
        with h5open(infos[0], "r") as h:
            mw.append(h.attrs["mag"])

    # plot
    x = []
    x.append(min(mw))
    x.append(max(mw))
    if x[0] < 6 < x[-1]:
        x.insert(1, 6)
    if x[0] < 8 < x[-1]:
        x.insert(-1, 8)

    fig_init()
    plt.plot(
        x, np.maximum(np.minimum((np.array(x) * 20 - 110), 50), 10), label="Target"
    )
    plt.plot(mw, nrup, label="Actual", marker="x", linestyle="None")

    plt.legend(loc="best")
    plt.title("Magnitude vs Number of Realizations")
    plt.ylabel("Number of Realizations")
    plt.xlabel("Magnitude")
    plt.savefig(os.path.join(out_dir, "mag_vs_nrup"))
    plt.close()


def plot_dbottom(srf_dirs, nhm_file, out_dir):
    srf_depths = []
    fault_names = []
    nhm_depths = []
    for d in srf_dirs:
        info = glob(os.path.join(d, "*.info"))[0]
        with h5open(info, "r") as i:
            srf_depths.append(np.average(i.attrs["dbottom"]))
        fault_names.append(os.path.basename(os.path.dirname(d)))
    # grab depths from nhm
    nhm_depths = np.zeros(len(fault_names))
    with open(nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
        n_i = 0
        while min(nhm_depths) <= 0:
            name = nhm[n_i].strip()
            if name in fault_names:
                nhm_depths[fault_names.index(name)] = float(nhm[n_i + 6].split()[0])
            elif name == "":
                print("Could not find all SRFs in NHM.")
                return
            n_i += 13 + int(nhm[n_i + 11])

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

    plt.legend(loc="best")
    plt.title("Seismogenic Depth from NHM and SRF")
    plt.ylabel("SRF Depth (km)")
    plt.xlabel("NHM Depth (km)")
    plt.savefig(os.path.join(out_dir, "dbottom"))
    plt.close()


def nhm2corners(nhm_file, names):
    corners_sub = []
    corners_shal = []
    ftypes = {}
    ex_names = []
    with open(nhm_file, "r") as n:
        for _ in range(15):
            n.readline()
        nhm = n.readlines()
    n_i = 0
    while n_i < len(nhm):
        n = nhm[n_i].strip()
        n_pt = int(nhm[n_i + 11])
        ftype = 2 + int(not nhm[n_i + 1].startswith("SUBDUCTION"))
        dbottom = float(nhm[n_i + 6].split()[0])
        dtop = float(nhm[n_i + 7].split()[0])
        if n not in names:
            ex_names.append(n)
            pwid = (dbottom - dtop + 3 * (dbottom >= 12)) / math.tan(
                math.radians(float(nhm[n_i + 3].split()[0]))
            )
            ftype -= 2
            corners = np.zeros((n_pt - 1, 4, 2))
            trace = nhm[n_i + 12 : n_i + 12 + n_pt]
            trace = [map(float, pair) for pair in map(str.split, trace)]
            dip_dir = float(nhm[n_i + 4])
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
        ftypes[n] = ftype
        n_i += 13 + n_pt

    return ftypes, corners_shal, corners_sub, ex_names


def plot_srf_nhm(srf_dirs, nhm_file, out_dir):
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
    faults = map(os.path.basename, map(os.path.dirname, srf_dirs))
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
    for i, d in enumerate(srf_dirs):
        name = faults[i]
        f = glob(os.path.join(d, "*.info"))[0]
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


def plot_hypo_dist(srf_dirs, out_dir):
    if os.path.isdir(out_dir):
        rmtree(out_dir)
    os.makedirs(out_dir)

    # normal distribution
    n_x = np.linspace(0, 1, 100)
    n_y = norm.cdf(n_x, 0.5, 0.25)
    n_call = lambda x: norm.cdf(x, 0.5, 0.25)
    # weibull distribution
    w_x = np.random.weibull(3.353, size=100000) * 0.612
    w_x.sort()
    w_y = np.arange(w_x.size) / (w_x.size - 1.0)
    w_call = lambda x: w_y[[np.argmin(np.abs(v - w_x)) for v in x]]

    # kstest
    names = []
    mw = []
    p_s = []
    p_d = []
    for d in srf_dirs:
        n = os.path.basename(os.path.dirname(d))
        names.append(n)
        r = glob(os.path.join(d, "*.info"))
        shypo = np.zeros(len(r) + 2)
        dhypo = np.zeros(len(r) + 2)
        shypo[-1] = 1
        dhypo[-1] = 1
        for i, f in enumerate(r):
            with h5open(f, "r") as h:
                shypo[i + 1] = h.attrs["shypo"][0] / sum(h.attrs["length"])
                dhypo[i + 1] = h.attrs["dhypo"][0] / h.attrs["width"][0]
                if not i:
                    mw.append(h.attrs["mag"])
        shypo.sort()
        dhypo.sort()
        y = np.arange(shypo.size) / (shypo.size - 1.0)

        # kstest
        p_s.append(kstest(shypo, n_call).pvalue)
        p_d.append(kstest(dhypo, w_call).pvalue)

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

        plt.legend(loc="best")
        plt.title(
            "%s Hypocentre Location Distribution (%d realizations)" % (n, shypo.size)
        )
        plt.ylabel("CDF")
        plt.xlabel("along fault")
        plt.gca().set_ylim([0, 1])
        plt.gca().set_xlim([0, 1])
        plt.savefig(os.path.join(out_dir, n))
        plt.close()

    # overall plot
    fig_init()
    plt.plot([min(mw), max(mw)], [0.05, 0.05])
    plt.plot(mw, p_s, label="p (along strike)", marker="x", linestyle="None")
    plt.plot(mw, p_d, label="p (along dip)", marker="x", linestyle="None")
    plt.legend(loc="best")
    plt.title("Hypocentre Location Distribution (%d simulations)" % (len(mw)))
    plt.ylabel("p")
    plt.xlabel("magnitude")
    plt.gca().set_yscale("log")
    # plt.gca().set_ylim([0, 1])
    plt.savefig(os.path.join(out_dir, "..", "hypo_dist"))
    plt.close()

    # overall histogram strike
    fig_init()
    plt.hist(p_s, bins=min(len(p_s), 50), label="along strike")
    plt.legend(loc="best")
    plt.title(
        "Hypocentre along Strike Distribution Probability of Fit (%d simulations)"
        % (len(mw))
    )
    plt.ylabel("n")
    plt.xlabel("probability of fitting normal distribution")
    plt.savefig(os.path.join(out_dir, "..", "hypo_dist_strike"))
    plt.close()
    # overall histogram dip
    fig_init()
    plt.hist(p_d, bins=min(len(p_s), 50), label="along dip")
    plt.legend(loc="best")
    plt.title(
        "Hypocentre along Dip Distribution Probability of Fit (%d simulations)"
        % (len(mw))
    )
    plt.ylabel("n")
    plt.xlabel("probability of fitting weibull distribution")
    plt.savefig(os.path.join(out_dir, "..", "hypo_dist_dip"))
    plt.close()

    # as text
    with open(os.path.join(out_dir, "..", "hypo_dist.txt"), "w") as d:
        d.write("name strike-p dip-p mw\n")
        for i in range(len(names)):
            d.write("%s %s %s %s\n" % (names[i], p_s[i], p_d[i], mw[i]))

###
### COMBINED PLOTS
###

def test_selection(selection_file, names, versus, out_dir):
    with open(os.path.join(out_dir, 'not_found_%s.txt' % (versus)), 'w') as o:
        with open(selection_file, "r") as s:
            for line in s:
                n = line.split()[0]
                if n not in names:
                    o.write("%s fault not found" % (n))

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


def plot_duration_mag(jsons, out_dir):
    durations = []
    magnitudes = []
    for j in jsons:
        with open(j, "r") as jo:
            d = json.load(jo)
        magnitudes.append(d["mag"])
        durations.append(d["sim_duration"])

    fig_init()
    plt.plot(
        magnitudes, durations, label="Velocity Models", marker="x", linestyle="None"
    )

    plt.legend(loc="best")
    plt.title("Magnitude vs Simulation Duration")
    plt.ylabel("Simulation Duration (s)")
    plt.xlabel("Magnitude")
    plt.savefig(os.path.join(out_dir, "mag_vs_duration"))
    plt.close()


def plot_size_mag(vm_dirs, out_dir):
    sizes = []
    magnitudes = []
    for d in vm_dirs:
        sizes.append(os.stat(os.path.join(d, "vs3dfile.s")).st_size / 1000000.0)
        with open(os.path.join(d, "params_vel.json"), "r") as jo:
            magnitudes.append(json.load(jo)["mag"])

    fig_init()
    plt.plot(
        magnitudes,
        sizes,
        label="Velocity Model Components",
        marker="x",
        linestyle="None",
    )

    plt.legend(loc="best")
    plt.title("Magnitude vs VM File Sizes")
    plt.ylabel("File Sizes (megabytes)")
    plt.xlabel("Magnitude")
    plt.savefig(os.path.join(out_dir, "mag_vs_size"))
    plt.close()


def plot_vm_params(jsons, out_dir):
    params = {"hh": [], "flo": []}
    for j in jsons:
        with open(j, "r") as jo:
            d = json.load(jo)
            params["hh"].append(d["hh"])
            params["flo"].append(d["flo"])
    x = np.arange(len(jsons))

    for parameter in params:
        fig_init()
        plt.plot(x, params[parameter], label=parameter, marker="x", linestyle="None")
        plt.legend(loc="best")
        plt.title("VM Parameter (%s)" % (parameter))
        plt.ylabel(parameter)
        plt.gca().get_xaxis().set_visible(False)
        plt.savefig(os.path.join(out_dir, "vm_%s" % (parameter)))
        plt.close()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("--srf_dir", help="NHM SRF directory to test")
    parser.add_argument("--vm_dir", help="NHM VM directory to test")
    parser.add_argument(
        "--nhm-file",
        help="NHM fault list file",
        default=os.path.join(os.path.dirname(__file__), "NZ_FLTmodel_2010.txt"),
    )
    parser.add_argument(
        "--selection-file", help="NHM selection file containing wanted faults"
    )
    parser.add_argument("--out-dir", help="Folder to place outputs in", default="./cstests")
    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    # run SRF plots and checks
    if args.srf_dir is not None:
        srf_dirs = glob(os.path.join(os.path.abspath(args.srf_dir), "*", "Srf"))
        srf_faults = map(os.path.basename, map(os.path.dirname, srf_dirs))
        info_files = glob(os.path.join(os.path.abspath(args.srf_dir), "*", "Srf", "*.info"))
        names = map(os.path.basename, map(os.path.dirname, srf_dirs))

        plot_area_mag(info_files, args.out_dir)
        plot_area_nhm(srf_dirs, args.nhm_file, args.out_dir)
        plot_mag_nhm(srf_dirs, args.nhm_file, args.out_dir)
        plot_mag_nrup(srf_dirs, args.out_dir)
        plot_dbottom(srf_dirs, args.nhm_file, args.out_dir)
        plot_srf_nhm(srf_dirs, args.nhm_file, args.out_dir)
        plot_hypo_dist(srf_dirs, os.path.join(args.out_dir, "hypo_dist"))
        if args.selection_file is not None:
            test_selection(args.selection_file, names, "SRF", args.out_dir)

    # run VM plots and checks
    if args.vm_dir is not None:
        vm_json = glob(os.path.join(args.vm_dir, "*", "params_vel.json"))
        vm_dirs = map(os.path.dirname, vm_json)
        names = map(os.path.basename, vm_dirs)

        plot_vm(vm_dirs, args.out_dir)
        plot_duration_mag(vm_json, args.out_dir)
        plot_size_mag(vm_dirs, args.out_dir)
        plot_vm_params(vm_json, args.out_dir)
        if args.selection_file is not None:
            test_selection(args.selection_file, names, "VM", args.out_dir)

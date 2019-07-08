#!/usr/bin/env python2
"""
REQUIREMENTS:
PATH contains plot_srf_square.py plot_srf_map.py if using -p (plot results)
PYTHONPATH contains createSRF and qcore
"""
import math
import os
import random
from subprocess import call
from time import time

import numpy as np
from scipy.stats import truncnorm

from createSRF import leonard, skarlatoudis, CreateSRF_multi, gen_meta
from qcore import geo, simulation_structure

SCRIPT_LOCATION = os.path.dirname(os.path.abspath(__file__))

MAGNITUDE_ROUNDING_THRESHOLD = 7.5


def get_seed():
    """Creates a seed for genslip.
    This number is the maximum value of a 32 bit int as that is how genslip stores the value"""
    return random.randint(2**32-1)


def rand_shyp_dhyp(length=1.0, width=1.0):
    # truncated normal distribution
    shyp_mu = 0.5
    shyp_sigma = 0.25
    lower_std_dev_limit = (0 - shyp_mu) / shyp_sigma
    upper_std_dev_limit = (1 - shyp_mu) / shyp_sigma

    shyp = truncnorm(
        lower_std_dev_limit, upper_std_dev_limit, loc=shyp_mu, scale=shyp_sigma
    ).rvs()

    # truncated weibull distribution
    dhyp_scale = 0.612
    dhyp_shape = 3.353

    # within range 0 -> 1 (exclusive)
    dhyp = np.random.weibull(dhyp_shape) * dhyp_scale
    while dhyp >= 1:
        dhyp = np.random.weibull(dhyp_shape) * dhyp_scale
    return shyp * length, dhyp * width


def nhm2info(nhm_file, nhm_skip, fault_name):
    dbi = nhm_skip
    with open(nhm_file, "r") as dbr:
        db = list(map(str.strip, dbr.readlines()))
    dbl = len(db)

    # Read through the NHM file to find the details for our fault
    while dbi < dbl and db[dbi] != fault_name:
        dbi += 13 + int(db[dbi + 11])

    # If we didn't find it return false
    if db[dbi] != fault_name:
        print("Didn't find fault in nhm")
        return False

    # Load data from NHM
    n_pt = int(db[dbi + 11])
    tect_type, fault_type = db[dbi + 1].split()
    # load trace
    traces = db[dbi + 12 : dbi + 12 + n_pt]
    pts = [list(map(float, ll.split())) for ll in traces]
    # clean points (remove duplicates)
    for i in range(n_pt - 2, -1, -1):
        if geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) < 0.1:
            del pts[i + 1]
            n_pt -= 1
    # total planes
    n_plane = n_pt - 1

    # other values from nhm file
    mag = [float(db[dbi + 10].split()[0])]
    rake = [[float(db[dbi + 5])] * n_plane]
    dip = [[float(db[dbi + 3].split()[0])] * n_plane]
    dtop = [[float(db[dbi + 7].split()[0])] * n_plane]
    dbottom = float(db[dbi + 6].split()[0])

    # strike orientation is based on dip direction
    dip_dir = float(db[dbi + 4])
    strike_avg = (dip_dir - 90) % 360

    return (
        n_pt,
        pts,
        mag,
        strike_avg,
        n_plane,
        dbottom,
        dtop,
        tect_type,
        rake,
        dip,
        dip_dir,
    )


def nhm2srf_gen_params(fault_name, rel_number, out_dir, nhm_skip, nhm_file):
    """Generates a dictionary that can be given to CreateSRF_multi with all required parameters"""
    # adjust outputs
    out_gmt = os.path.join(out_dir, "fault_traces.gmt")
    out_log = os.path.join(out_dir, "logfile.txt")
    # prepare file system
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(out_gmt):
        os.remove(out_gmt)
    if not os.path.exists(out_log):
        with open(out_log, "w") as log:
            log.write("filename\tshyp\tdhyp\tseed\n")

    # load db
    n_pt, pts, mag, strike_avg, n_plane, dbottom, dtop, tect_type, rake, dip, dip_dir = nhm2info(
        nhm_file, nhm_skip, fault_name
    )
    # Finished loading from NHM

    # midpoints, lengths and strikes of all faults
    mids = []
    lengths = []
    strikes = []
    stk_norm = 0
    stk_rev = 0

    for plane in range(n_pt - 1):
        mids.append(
            geo.ll_mid(
                pts[plane][0], pts[plane][1], pts[plane + 1][0], pts[plane + 1][1]
            )
        )
        dist = geo.ll_dist(
            pts[plane][0], pts[plane][1], pts[plane + 1][0], pts[plane + 1][1]
        )
        lengths.append(round_subfault_size(dist, mag))
        bearing = geo.ll_bearing(
            mids[plane][0], mids[plane][1], pts[plane + 1][0], pts[plane + 1][1]
        )
        if abs((bearing - strike_avg + 180) % 360 - 180) < 90:
            strikes.append(bearing)
            stk_norm += 1
        else:
            strikes.append((bearing + 180) % 360)
            stk_rev += 1

    # Check direction of the fault
    if stk_norm and stk_rev:
        print(
            "WARNING: FAULT GOES BACK IN REVERSE OF ORIGINAL DIRECTION: {}".format(
                fault_name
            )
        )
        return False

    # assuming no reverse angles and ordering in one direction
    if stk_rev:
        mids = mids[::-1]
        lengths = lengths[::-1]
        strikes = strikes[::-1]
    trace_length = sum(lengths)

    # fixed values
    dt = 0.025
    cases = ["combined"]
    nseg = [n_plane]
    seg_delay = [0]
    mom = [-1]
    rvfac_seg = ["-1"]
    gwid = ["-1"]
    rup_delay = [0]

    flen = [lengths]
    # subfault density
    dlen = [[0.1] * n_plane]
    # Karim: add 3km to depth if bottom >= 12km
    if dbottom >= 12:
        raw_fwid = [
            [(dbottom - dtop[0][0] + 3) / math.sin(math.radians(dip[0][0]))] * n_plane
        ]
    else:
        raw_fwid = [
            [(dbottom - dtop[0][0]) / math.sin(math.radians(dip[0][0]))] * n_plane
        ]

    fwid = [[round_subfault_size(f, mag) for f in segment] for segment in raw_fwid]

    # Scale magnitude based on type of interface
    if tect_type == "SUBDUCTION_INTERFACE":
        mag = [skarlatoudis(fwid[0][0] * trace_length)]
    else:
        mag = [leonard(rake[0][0], fwid[0][0] * trace_length)]

    stk = [strikes]
    elon = [[ll[0] for ll in mids]]
    elat = [[ll[1] for ll in mids]]

    shyp_shift, dhyp_shift = rand_shyp_dhyp(trace_length, fwid[0][0])

    shypo = [[shyp_shift - (lengths[0] / 2.0)]]
    dhypo = [[dhyp_shift] * n_plane]
    prefix = os.path.join(
        out_dir,
        simulation_structure.get_srf_location(
            simulation_structure.get_realisation_name(fault_name, rel_number)
        ),
    )[:-4]

    seed = get_seed()

    return {
        "nseg": nseg,
        "seg_delay": seg_delay,
        "mag0": mag,
        "mom0": mom,
        "rvfac_seg": rvfac_seg,
        "gwid": gwid,
        "rup_delay": rup_delay,
        "flen": flen,
        "dlen": dlen,
        "fwid": fwid,
        "dwid": dlen,
        "dtop": dtop,
        "stk": stk,
        "rak": rake,
        "dip": dip,
        "elon": elon,
        "elat": elat,
        "shypo": shypo,
        "dhypo": dhypo,
        "dt": dt,
        "seed": seed,
        "prefix0": prefix,
        "cases": cases,
        "dip_dir": dip_dir,
        "silent": False,
        "tect_type": tect_type,
    }


def generate_srf(out_dir, fault_name, create_srf_params, plot):
    """Takes the srf params dictionary and uses it to generate """
    # store parameters
    out_log = os.path.join(out_dir, "logfile.txt")
    with open(out_log, "a") as log:
        log.write(
            "{}.srf\t{}\t{}\t{}\n".format(
                create_srf_params["prefix0"],
                create_srf_params["shypo"][0][0],
                create_srf_params["dhypo"][0][0],
                create_srf_params["seed"],
            )
        )

    # store fault traces
    out_gmt = os.path.join(out_dir, "fault_traces.gmt")
    with open(out_gmt, "a") as traces:
        traces.write("> {}\n{}\n".format(fault_name, "\n".join(traces)))

    t0 = time()
    print("creating SRF: {}".format(fault_name))
    # Splat the parameters dictionary into createSRF_multi
    CreateSRF_multi(**create_srf_params)
    print("created SRF: {} {:.2f}s)".format(fault_name, time() - t0))

    if plot:
        t0 = time()
        print("plotting SRF: {}".format(fault_name))
        call(["plot_srf_square.py", "{}.srf".format(create_srf_params["prefix0"])])
        call(["plot_srf_map.py", "{}.srf".format(create_srf_params["prefix0"])])
        print("plotted SRF: {} ({:.2f}s)".format(fault_name, time() - t0))

    srf_file = os.path.join(
        out_dir, fault_name, "Srf", "{}_REL01.srf".format(fault_name)
    )

    gen_meta(
        srf_file,
        4,
        create_srf_params["mag0"],
        create_srf_params["stk"],
        create_srf_params["rake"],
        create_srf_params["dip"],
        create_srf_params["dt"],
        vm="{}/lp_generic1d-gp01_v1.vmod".format(
            os.path.dirname(os.path.abspath(__file__))
        ),
        dip_dir=create_srf_params["dip_dir"],
        shypo=[
            s[0]
            + 0.5 * create_srf_params["flen"][len(create_srf_params["cases"]) - 1][0]
            for s in create_srf_params["shypo"]
        ],
        dhypo=[d[0] for d in create_srf_params["dhypo"]],
        tect_type=create_srf_params["tect_type"],
        file_name=os.path.join(out_dir, fault_name, fault_name),
    )

    return True


def round_subfault_size(dist, mag):
    if mag[0] > MAGNITUDE_ROUNDING_THRESHOLD:
        return round(dist * 2) / 2
    else:
        return round(dist * 10) / 10


def main():
    from argparse import ArgumentParser

    # parameters
    parser = ArgumentParser()
    parser.add_argument("fault", help="The name of the fault to be generated")
    parser.add_argument(
        "realisation",
        nargs="?",
        default=1,
        type=int,
        help="The number of the realisation to be generated",
    )
    parser.add_argument(
        "-o", "--out_dir", help="directory to place outputs", default="Sources"
    )
    parser.add_argument(
        "--nhm_file",
        help="NHM file location",
        default=os.path.join(SCRIPT_LOCATION, "NZ_FLTmodel_2010.txt"),
    )  # This file is a symbolic link: To avoid need to edit hard-coded filename
    parser.add_argument(
        "--nhm_skip", help="NHM header lines to skip", type=int, default=15
    )
    parser.add_argument(
        "--dhypo", default=[0.6], type=float, help="depth of hypocentre"
    )
    parser.add_argument("-p", "--plot", help="plot results", action="store_true")
    args = parser.parse_args()
    out_dir = os.path.abspath(args.out_dir)

    fault_name = args.fault

    # load wanted fault information
    create_srf_params = nhm2srf_gen_params(
        fault_name, args.realisation, out_dir, args.nhm_skip, args.nhm_file
    )
    generate_srf(out_dir, fault_name, create_srf_params, args.plot)


if __name__ == "__main__":
    main()

#!/usr/bin/env python2
"""
REQUIREMENTS:
PATH contains plot_srf_square.py plot_srf_map.py if using -p (plot results)
PYTHONPATH contains createSRF and qcore
"""

import math
import os
from logging import Logger
from multiprocessing.pool import Pool
from subprocess import call
import sys
from time import time
import logging

import numpy as np

from qcore import geo, simulation_structure, qclogging

from createSRF import leonard, skarlatoudis, CreateSRF_multi, gen_meta

MAGNITUDE_ROUNDING_THRESHOLD = 7.5


def rand_shyp_dhyp(length=1.0, width=1.0):
    # normal distribution
    shyp_mu = 0.5
    shyp_sigma = 0.25
    # weibel distribution
    dhyp_scale = 0.612
    dhyp_shape = 3.353

    shyp = -1
    dhyp = 2

    # within range 0 -> 1 (exclusive)
    while shyp <= 0 or shyp >= 1:
        shyp = np.random.normal(shyp_mu, shyp_sigma)
    while dhyp >= 1:
        dhyp = np.random.weibull(dhyp_shape) * dhyp_scale
    return shyp * length, dhyp * width


###
### PREPARE TASKS
###
# most work is trivial, only need to run CreateSRF_multi parallel
def load_msgs(args, fault_names, faults, logger: Logger = qclogging.get_basic_logger()):
    # adjust outputs
    out_gmt = os.path.join(args.out_dir, "fault_traces.gmt")
    out_log = os.path.join(args.out_dir, "logfile.txt")
    # prepare file system
    if not os.path.isdir(args.out_dir):
        logger.debug("Making out_dir: {}".format(args.out_dir))
        os.makedirs(args.out_dir)
        qclogging.remove_buffer_handler(logger)
    if os.path.exists(out_gmt):
        logger.debug("Removing GMT directory: {}".format(out_gmt))
        os.remove(out_gmt)
    logger.debug("Using {} as the srfgen variable log".format(out_log))
    with open(out_log, "w") as log:
        log.write("filename\tshyp\tdhyp\tseed\n")

    # load db
    dbi = args.nhm_skip
    with open(args.nhm_file, "r") as dbr:
        db = list(map(str.strip, dbr.readlines()))
    dbl = len(db)

    # compile messages (parameters for CreateSRF_multi)
    if fault_names is None:
        return []

    msgs = []
    while dbi < dbl:
        name = db[dbi]
        tect_type, fault_type = db[dbi + 1].split()
        # points in fault trace
        n_pt = int(db[dbi + 11])
        skip = 13 + n_pt
        # skip if not wanted
        if name not in fault_names:
            logging.log(
                logging.DEBUG // 2,
                "{} not in the list of fault names, continuing".format(name),
            )
            dbi += skip
            continue
        logger.debug("Found {}, procceding to add to messages".format(name))

        # load trace
        pts = [list(map(float, ll.split())) for ll in db[dbi + 12 : dbi + 12 + n_pt]]
        # clean points (remove duplicates)
        for i in range(n_pt - 2, -1, -1):
            if geo.ll_dist(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]) < 0.1:
                del pts[i + 1]
                n_pt -= 1
        # total planes
        n_plane = n_pt - 1
        # strike orientation is based on dip direction
        dip_dir = float(db[dbi + 4])
        strike_avg = (dip_dir - 90) % 360
        # midpoints, lengths and strikes of all faults
        mids = []
        lengths = []
        strikes = []
        stk_norm = 0
        stk_rev = 0

        # other values from nhm file
        mag = [float(db[dbi + 10].split()[0])]
        rake = [[float(db[dbi + 5])] * n_plane]
        dip = [[float(db[dbi + 3].split()[0])] * n_plane]
        dtop = [[float(db[dbi + 7].split()[0])] * n_plane]

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
        if stk_norm and stk_rev:
            logger.critical(
                "WARNING: FAULT GOES BACK IN REVERSE OF ORIGINAL DIRECTION: {}".format(
                    name
                )
            )
        # assuming no reverse angles and ordering in one direction
        if stk_rev:
            mids = mids[::-1]
            lengths = lengths[::-1]
            strikes = strikes[::-1]
        trace_length = sum(lengths)

        # wanted parameters to override
        t_hypo = "n"
        n_hypo = args.nhypo
        n_slip = args.nslip
        dhypos = args.dhypo
        if faults is not None:
            fault = faults[fault_names.index(name)]
            if len(fault) >= 2:
                # given as hypocentre every x km
                if fault[1][-1] == "k":
                    t_hypo = "k"
                    hyp_step = float(fault[1][:-1])
                    n_hypo = 1 + int(trace_length // hyp_step)
                    # 0th hypocentre position
                    if n_hypo == 1:
                        z_hypo = trace_length / 2.0
                    else:
                        z_hypo = (trace_length % hyp_step) / 2.0
                    logger.debug(
                        "Fault found in fault selection file, making {} realisations with equally spaced hypocenters".format(
                            hyp_step
                        )
                    )

                elif fault[1][-1] == "n":
                    t_hypo = "n"
                    n_hypo = int(fault[1][:-1])
                    hyp_step = trace_length / (n_hypo * 2.0)

                # given as number of randomly placed hypocentres
                elif fault[1][-1] == "r":
                    t_hypo = "r"
                    n_hypo = int(fault[1][:-1])
                    logger.debug(
                        "Fault found in fault selection file. Making {} realisations with randomised parameters".format(
                            n_hypo
                        )
                    )

                # given as number of hypocentres
                else:
                    n_hypo = int(fault[1])
                    logger.debug(
                        "Fault found in fault selection file. Making {} realisations".format(
                            n_hypo
                        )
                    )

            if len(fault) >= 3:
                n_slip = int(fault[2])
                logger.debug("n_slip set to {}".format(n_slip))
            if len(fault) >= 4:
                dhypos = list(map(float, fault[3].split(",")))
                logger.debug("dhypos set to {}".format(dhypos))
        if t_hypo == "n":
            hyp_step = trace_length / (n_hypo * 2.0)
            logger.debug("hyp_step set to {}".format(hyp_step))
        seed = args.seed

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
        if float(db[dbi + 6].split()[0]) >= 12:
            raw_fwid = [
                [
                    (float(db[dbi + 6].split()[0]) - dtop[0][0] + 3)
                    / math.sin(math.radians(dip[0][0]))
                ]
                * n_plane
            ]
        else:
            raw_fwid = [
                [
                    (float(db[dbi + 6].split()[0]) - dtop[0][0])
                    / math.sin(math.radians(dip[0][0]))
                ]
                * n_plane
            ]
        fwid = []
        for segment in raw_fwid:
            fwid.append([round_subfault_size(f, mag) for f in segment])

        if tect_type == "SUBDUCTION_INTERFACE":
            logger.debug(
                "Subduction interface tectonic type, using Skarlatoudis formula to determine magnitude"
            )
            mag = [skarlatoudis(fwid[0][0] * trace_length)]
        else:
            logger.debug(
                "Non subduction interface tectonic type, using Leonard formula to determine magnitude"
            )
            mag = [leonard(rake[0][0], fwid[0][0] * trace_length)]
        dwid = dlen
        stk = [strikes]
        elon = [[ll[0] for ll in mids]]
        elat = [[ll[1] for ll in mids]]

        for n_shyp in range(n_hypo):
            # hypocentre position from far left edge
            if t_hypo == "n":
                shyp_shift = hyp_step * (1 + 2 * n_shyp)
            elif t_hypo == "k":
                shyp_shift = z_hypo + hyp_step * n_shyp
            elif t_hypo == "r":
                shyp_shift, dhyp_shift = rand_shyp_dhyp(trace_length, fwid[0][0])
            logger.debug("shyp_shift set to: {}".format(shyp_shift))
            # NOTE: this shypo is relative to the first combined fault
            # if not adjusted later, must be relative to full length
            shypo = [[shyp_shift - (lengths[0] / 2.0)]]
            for _ in range(n_slip):
                for i, d in enumerate(dhypos):
                    seed += args.seed_inc
                    if t_hypo == "r":
                        dhypo = [[dhyp_shift] * n_plane]
                    else:
                        dhypo = [[fwid[0][0] * d] * n_plane]
                    prefix = os.path.join(
                        args.out_dir,
                        simulation_structure.get_srf_location(
                            simulation_structure.get_realisation_name(name, n_shyp + 1)
                        ),
                    )[:-4]
                    # create SRF from description
                    msgs.append(
                        {
                            "nseg": nseg,
                            "seg_delay": seg_delay,
                            "mag": mag,
                            "mom": mom,
                            "rvfac_seg": rvfac_seg,
                            "gwid": gwid,
                            "rup_delay": rup_delay,
                            "flen": flen,
                            "dlen": dlen,
                            "fwid": fwid,
                            "dwid": dwid,
                            "dtop": dtop,
                            "stk": stk,
                            "rake": rake,
                            "dip": dip,
                            "elon": elon,
                            "elat": elat,
                            "shypo": shypo,
                            "dhypo": dhypo,
                            "dt": dt,
                            "seed": seed,
                            "prefix": prefix,
                            "cases": cases,
                            "dip_dir": dip_dir,
                            "stoch": os.path.join(args.out_dir, name, "Stoch"),
                            "name": name,
                            "tect_type": tect_type,
                            "plot": args.plot,
                            "genslip_version": "3.3",
                        }
                    )
                    if tect_type == "SUBDUCTION_INTERFACE":
                        msgs[-1].update({"genslip_version": "5.4.2"})
                    # store parameters
                    with open(out_log, "a") as log:
                        log.write(
                            "%s.srf\t%s\t%s\t%s\n"
                            % (prefix, shypo[0][0], dhypo[0][0], seed)
                        )
        # Initialise the logger so each process can get the local copy
        f_logger = qclogging.get_realisation_logger(logger, name)
        f_logger.propagate = False

        # store fault traces
        with open(out_gmt, "a") as traces:
            traces.write(
                "> %s\n%s\n" % (name, "\n".join(db[dbi + 12 : dbi + 12 + n_pt]))
            )

        # move to next fault definition
        dbi += skip
    logger.debug("{} messages collected".format(len(msgs)))
    return msgs


def round_subfault_size(dist, mag):
    if mag[0] > MAGNITUDE_ROUNDING_THRESHOLD:
        return round(dist * 2) / 2
    else:
        return round(dist * 10) / 10


def run_create_srf(fault):
    logger = qclogging.get_logger(fault["name"])
    t0 = time()
    #    sys.stdout = open(str(os.getpid())+".out","w")
    logger.info("Creating SRF: {}".format(fault["name"]))
    # all of the work, rest of the script is complete under 1 second
    CreateSRF_multi(
        fault["nseg"],
        fault["seg_delay"],
        fault["mag"],
        fault["mom"],
        fault["rvfac_seg"],
        fault["gwid"],
        fault["rup_delay"],
        fault["flen"],
        fault["dlen"],
        fault["fwid"],
        fault["dwid"],
        fault["dtop"],
        fault["stk"],
        fault["rake"],
        fault["dip"],
        fault["elon"],
        fault["elat"],
        fault["shypo"],
        fault["dhypo"],
        fault["dt"],
        fault["seed"],
        fault["prefix"],
        fault["cases"],
        dip_dir=fault["dip_dir"],
        stoch=fault["stoch"],
        tect_type=fault["tect_type"],
        silent=True,
        genslip_version=fault["genslip_version"],
        logger=logger,
    )
    logger.debug("created SRF: {} ({:.2f}s)".format(fault["name"], time() - t0))
    if fault["plot"]:
        t0 = time()
        logger.info("plotting SRF: {}".format(fault["name"]))
        call(["plot_srf_square.py", "{}.srf".format(fault["prefix"])])
        call(["plot_srf_map.py", "{}.srf".format(fault["prefix"])])
        logger.info("plotted SRF: {} (:.2fs)".format(fault["name"], time() - t0))
    srf_file = os.path.join(
        args.out_dir, fault["name"], "Srf", "{}_REL01.srf".format(fault["name"])
    )
    gen_meta(
        srf_file,
        4,
        fault["mag"],
        fault["stk"],
        fault["rake"],
        fault["dip"],
        fault["dt"],
        vm="%s/lp_generic1d-gp01_v1.vmod"
        % (os.path.dirname(os.path.abspath(__file__))),
        dip_dir=fault["dip_dir"],
        shypo=[
            s[0] + 0.5 * fault["flen"][len(fault["cases"]) - 1][0]
            for s in fault["shypo"]
        ],
        dhypo=[d[0] for d in fault["dhypo"]],
        tect_type=fault["tect_type"],
        file_name=os.path.join(args.out_dir, fault["name"], fault["name"]),
        logger=logger,
    )


#   sys.stdout.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    logger = qclogging.get_logger("nhm2srf")
    # parameters
    parser = ArgumentParser()
    arg = parser.add_argument
    arg("selection_file", help="fault selection file")
    parser.add_argument(
        "-o", "--out-dir", help="directory to place outputs", default="autosrf"
    )
    arg(
        "--nhm-file",
        help="NHM file location",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "NZ_FLTmodel_2010.txt"
        ),
    )  # This file is a symbolic link: To avoid need to edit hard-coded filename
    arg("--nhm-skip", help="NHM header lines to skip", type=int, default=15)
    arg("--seed", help="initial seed", type=int, default=1234)
    arg("--seed-inc", help="seed increment", type=int, default=10)
    arg("--nhypo", default=1, help="hypocentre number/spacing if not in selection_file")
    arg("--nslip", type=int, default=1, help="number of slips if not in selection_file")
    arg(
        "--dhypo",
        action="append",
        type=float,
        help="depth of hypocentre, repeat parameter as needed",
    )
    arg("-n", "--nproc", help="number of processes", type=int, default=1)
    arg("-p", "--plot", help="plot results", action="store_true")
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)
    if os.path.isdir(args.out_dir):
        qclogging.add_general_file_handler(logger, os.path.join("nhm2srf_log.txt"))
    else:
        qclogging.add_buffer_handler(logger, file_name=os.path.join("nhm2srf_log.txt"))

    # selection file or ALL
    if args.selection_file == "ALL":
        faults = None
        fault_names = None
        logger.debug("Creating realisations for all faults")
    elif not os.path.exists(args.selection_file):
        logger.error("Fault selecion file not found: {}".format(args.selection_file))
        sys.exit(1)
    else:
        with open(args.selection_file, "r") as select:
            faults = list(map(str.split, select.readlines()))
        fault_names = [f[0] for f in faults]
        logger.debug(
            "Creating realisations for the following faults: {}".format(
                ", ".join(fault_names)
            )
        )
    # default value
    if args.dhypo is None:
        args.dhypo = [0.6]
        logger.debug("dhypo not set, using {}".format(args.dhypo))

    # load wanted fault information
    msg_list = load_msgs(args, fault_names, faults, logger=logger)
    if len(msg_list) == 0:
        logger.info("No matches found in the NHM file, exiting.")
        sys.exit(1)

    # distribute work
    p = Pool(args.nproc)
    p.map(run_create_srf, msg_list)

#!/usr/bin/env python3

from argparse import ArgumentParser
import os
from logging import Logger
from subprocess import call

import numpy as np
import pandas as pd
import yaml

from qcore.pool_wrapper import PoolWrapper
from qcore import qclogging

from SrfGen.createSRF import CreateSRF_ps, gen_meta
from SrfGen.createSourceRealisation import create_ps_realisation


def run_create_srf(
    args, t, vs, rho, n_sims, logger: Logger = qclogging.get_basic_logger()
):
    mom = -1
    try:
        pid = t.pid.decode()  # t.pid is originally numpy.bytes_
    except AttributeError:
        pid = t.pid
    logger.debug(
        "Entered run_create_srf. pid is {}. Switching to fault logger".format(pid)
    )
    fault_logger = qclogging.get_realisation_logger(logger, pid)

    prefix = os.path.join(args.out_dir, pid, "Srf", pid)
    fault_logger.debug("Using prefix: {}".format(prefix))
    if args.uncertainty_file:
        fault_logger.debug(
            "Uncertainty file found. Generating perturbated realisations"
        )
        with open(args.uncertainty_file) as ao_f:
            add_opts = yaml.load(ao_f)
        create_ps_realisation(
            args.out_dir,
            pid,
            t.lat,
            t.lon,
            t.depth,
            t.mag,
            mom,
            t.strike,
            t.rake,
            t.dip,
            n_realisations=n_sims,
            additional_options=add_opts,
            dt=args.dt,
            vs=vs,
            rho=rho,
            silent=True,
            logger=fault_logger,
        )
        srf_file = os.path.join(args.out_dir, pid, "Srf", "{}_REL01.srf".format(pid))
        gen_meta(
            srf_file,
            1,
            t.mag,
            t.strike,
            t.rake,
            t.dip,
            0.005,
            lon=t.lon,
            lat=t.lat,
            vs=vs,
            rho=rho,
            centroid_depth=t.depth,
            file_name=os.path.join(args.out_dir, pid, pid),
            logger=fault_logger,
        )
    else:
        stoch = os.path.join(args.out_dir, pid, "Stoch")
        CreateSRF_ps(
            t.lat,
            t.lon,
            t.depth,
            t.mag,
            mom,
            t.strike,
            t.rake,
            t.dip,
            dt=args.dt,
            prefix=prefix,
            stoch=stoch,
            vs=vs,
            rho=rho,
            silent=True,
            logger=fault_logger,
        )
    if args.plot:
        fault_logger.debug("Plotting srf map")
        call(["plot_srf_map.py", "%s.srf" % (prefix)])


def run_create_srf_star(args_t_vs_rho):
    return run_create_srf(*args_t_vs_rho)


if __name__ == "__main__":
    # load parameters
    logger = qclogging.get_logger("gcmt2srf")
    parser = ArgumentParser()
    parser.add_argument("csv_file", help="path to CMT solutions CSV")
    parser.add_argument("velocity_model", help="path to 1D velocity model")
    parser.add_argument(
        "-o", "--out-dir", help="directory to place outputs", default="./autosrf"
    )
    parser.add_argument("--dt", help="SRF timestep", type=float, default=0.005)
    parser.add_argument("-n", "--nproc", type=int, default=1)
    parser.add_argument("-p", "--plot", action="store_true")
    uncertainty_group = parser.add_argument_group(
        "Uncertainty options",
        "Options to be used when SRFs are to be generated with uncertainty factors",
    )
    uncertainty_group.add_argument(
        "-r", "--realistion-uncertainty", dest="uncertainty_file", type=str, default=""
    )
    uncertainty_group.add_argument(
        "-c", "--cybershake-fault-file", dest="cs_file", type=str, default=""
    )
    args = parser.parse_args()
    # validate parameters
    qclogging.add_general_file_handler(
        logger, os.path.join(args.out_dir, "gcmt2srf_log.txt")
    )

    error_messages = []
    if not os.path.exists(args.csv_file):
        error_messages.append(
            "CMT solutions CSV file not found: {}".format(args.csv_file)
        )

    if not os.path.exists(args.velocity_model):
        error_messages.append(
            "Velocity model file not found: {}".format(args.velocity_model)
        )

    all_opts = bool(args.uncertainty_file) and bool(args.cs_file)
    any_opts = bool(args.uncertainty_file) or bool(args.cs_file)
    if any_opts and not all_opts:
        error_messages.append(
            "If realisation uncertainty or a cybershake fault file is given then the other must be given also."
        )

    if len(error_messages) > 0:
        for message in error_messages:
            logger.log(qclogging.NOPRINTERROR, message)
        parser.error("\n".join(error_messages))

    # load CMT CSV
    sources = np.rec.array(
        np.loadtxt(
            args.csv_file,
            dtype=[
                ("pid", "|S32"),
                ("lat", "f8"),
                ("lon", "f8"),
                ("depth", "f8"),
                ("mag", "f8"),
                ("strike", "f8"),
                ("dip", "f8"),
                ("rake", "f8"),
            ],
            delimiter=",",
            skiprows=1,
            usecols=(0, 2, 3, 13, 11, 4, 5, 6),
        )
    )
    logger.debug("Loaded {} sources".format(len(sources)))

    if all_opts:
        # Sort the array by pid to match the pandas array
        sources = np.sort(sources, order="pid")

        with open(args.cs_file) as cs_f:
            cs_file_pd = pd.read_csv(
                cs_f, sep=" ", header=None, names=["pid", "N_Rels"]
            )

        # Apply the identity function to the PublicID column and a filter function to the second column
        cs_file_pd = cs_file_pd.transform(
            {
                "pid": lambda x: x,
                "N_Rels": lambda x: int(x[:-1]) if x.endswith("r") else 0,
            }
        )

        # Remove all values that failed the filter function
        cs_file_pd = cs_file_pd[cs_file_pd.N_Rels != 0]

        n_sims = list(cs_file_pd.sort_values("pid")["N_Rels"])
    else:
        logger.debug(
            "No cybershake file given, creating one realisation for each source"
        )
        n_sims = np.ones(len(sources))

    # load velocity model (vs and rho)
    vmodel = np.rec.array(
        np.loadtxt(
            args.velocity_model,
            skiprows=1,
            usecols=(0, 2, 3),
            dtype=[("depth", "f8"), ("vs", "f8"), ("rho", "f8")],
        )
    )
    depth_bins = np.digitize(
        np.around(sources.depth, decimals=5),
        np.around(np.cumsum(vmodel.depth), decimals=5),
    )
    vs = vmodel.vs[depth_bins]
    rho = vmodel.rho[depth_bins]

    # distribute work
    p = PoolWrapper(args.nproc)
    p.map(
        run_create_srf_star, list(zip([args] * len(sources), sources, vs, rho, n_sims))
    )

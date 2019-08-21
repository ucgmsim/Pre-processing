#!/usr/bin/env python3

from argparse import ArgumentParser
import os
from logging import Logger
from multiprocessing.pool import Pool
from subprocess import call
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import yaml

from qcore import qclogging
from createSRF import (
    leonard,
    CreateSRF_ps,
    gen_meta,
    create_srf_psff,
    focal_mechanism_2_finite_fault,
)
from createSourceRealisation import create_ps_realisation, create_ps2ff_realisation

from createSRF import CreateSRF_ps, gen_meta
from createSourceRealisation import create_ps_realisation

def run_create_ps_srf(
    args, t, vs, rho, n_sims, logger_name: str = "Basic"
):
    logger = qclogging.get_logger(logger_name)
    mom = -1
    try:
        pid = t.pid.decode()  # t.pid is originally numpy.bytes_
    except AttributeError:
        pid = t.pid
    logger.debug(
        "Entered run_create_srf. pid is {}. Switching to fault logger".format(pid)
    )
    fault_logger = qclogging.get_realisation_logger(logger, pid)
    fault_logger.propagate = False

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


def run_create_psff_srf(args, t, n_sims):
    try:
        pid = t.pid.decode()  # t.pid is originally numpy.bytes_
    except AttributeError:
        pid = t.pid

    prefix = os.path.join(args.out_dir, pid, "Srf", pid)
    if args.uncertainty_file:
        with open(args.uncertainty_file) as ao_f:
            add_opts = yaml.load(ao_f)
        create_ps2ff_realisation(
            args.out_dir,
            pid,
            t.lat,
            t.lon,
            t.depth,
            t.mag,
            t.strike,
            t.rake,
            t.dip,
            n_realisations=n_sims,
            additional_options=add_opts,
            dt=args.dt,
            mwsr=args.mwsr,
            seed=args.seed,
            genslip_version=args.genslip_version,
            silent=True,
        )
        srf_file = os.path.join(args.out_dir, pid, "Srf", "{}_REL01.srf".format(pid))
        shypo, dhypo = focal_mechanism_2_finite_fault(
            t.lat, t.lon, t.depth, t.mag, t.strike, t.rake, t.dip, args.mwsr
        )[-2:]
        gen_meta(
            srf_file,
            2,
            t.mag,
            t.strike,
            t.rake,
            t.dip,
            args.dt,
            lon=t.lon,
            lat=t.lat,
            vm=args.velocity_model,
            shypo=shypo,
            dhypo=dhypo,
            file_name=os.path.join(args.out_dir, pid, pid),
            mwsr=args.mwsr,
        )
    else:
        stoch = os.path.join(args.out_dir, pid, "Stoch")
        create_srf_psff(
            t.lat,
            t.lon,
            t.mag,
            t.strike,
            t.rake,
            t.dip,
            t.depth,
            args.dt,
            prefix,
            seed=args.seed,
            stoch=stoch,
            genslip_version=args.genslip_version,
            mwsr=args.mwsr,
        )
    if args.plot:
        call(["plot_srf_map.py", "%s.srf" % (prefix)])


def main():
    # load parameters
    parser = ArgumentParser(add_help=False)
    parser.add_argument("csv_file", help="path to CMT solutions CSV")
    parser.add_argument("velocity_model", help="path to 1D velocity model")
    parser.add_argument(
        "type",
        help="The type of fault for the source being produced. 1 for point source, 2 for finite fault from point source",
        choices=[1, 2],
        type=int,
    )
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
    args, extra_args = parser.parse_known_args()
    if args.type == 1:
        pass
    elif args.type == 2:
        psff_parser = ArgumentParser(parents=[parser])
        psff_parser.add_argument(
            "-m",
            "--mwsr",
            help="Magnitude scaling relation. Only used for type 2 faults.",
        )
        psff_parser.add_argument(
            "-s" "--seed",
            help="The seed to be passed to the srf generation binary. Only used for type 2 faults.",
            default=None,
            type=int,
        )
        psff_parser.add_argument(
            "-g", "--genslip_version", type=str, default="3.3", options=["3.3", "5.4"]
        )
        args = psff_parser.parse_args()

    # validate parameters
    if not os.path.exists(args.csv_file):
        sys.exit("CMT solutions CSV file not found: %s" % (args.csv_file))
    if not os.path.exists(args.velocity_model):
        sys.exit("Velocity model file not found: %s" % (args.velocity_model))
    all_opts = bool(args.uncertainty_file) and bool(args.cs_file)
    any_opts = bool(args.uncertainty_file) or bool(args.cs_file)
    if any_opts and not all_opts:
        parser.error(
            "If realisation uncertainty or a cybershake fault file is given then the other must be given also."
        )
        exit(1)
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
    out_dir = os.path.abspath(args.out_dir)
    if os.path.isdir(out_dir):
        qclogging.add_general_file_handler(
            logger, os.path.join(args.out_dir, "gcmt2srf_log.txt")
        )
    else:
        qclogging.add_buffer_handler(logger, file_name=os.path.join(args.out_dir, "gcmt2srf_log.txt"))

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

    if args.type == 1:
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

        source_gen_function = run_create_ps_srf
        messages = list(zip([args] * len(sources), sources, vs, rho, n_sims))
    elif args.type == 2:
        source_gen_function = run_create_psff_srf
        messages = list(zip([args] * len(sources), sources, n_sims))
    else:
        source_gen_function = lambda: 0
        messages = []

    # distribute work
    p = Pool(args.nproc)
    p.starmap(source_gen_function, messages)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example (set out_dir with -o or --out-dir):
mpirun -n 3 srfinfo2vm.py "Srf/*.info"

HELP:
./srfinfo2vm -h
"""

from argparse import ArgumentParser
import csv
from distutils.spawn import find_executable
from glob import glob
from h5py import File as h5open
from logging import Logger
import math
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import platform
import subprocess
from shutil import rmtree, move, copyfile
from subprocess import Popen
import sys
from tempfile import mkdtemp
import yaml

from srfinfo2vm import create_vm,store_summary

from srf_generation.input_file_generation.realisation_to_srf import get_corners
from qcore import constants, geo, gmt, qclogging
from qcore.geo import R_EARTH
from qcore.utils import dump_yaml
from qcore.validate_vm import validate_vm

from empirical.util.classdef import GMM, Site, Fault
from empirical.util.empirical_factory import compute_gmm

from createSRF import leonard, mag2mom, mom2mag
from gen_coords import gen_coords

script_dir = os.path.dirname(os.path.abspath(__file__))
NZ_CENTRE_LINE = os.path.join(script_dir, "NHM/res/centre.txt")
NZ_LAND_OUTLINE = os.path.join(script_dir, "NHM/res/rough_land.txt")
NZVM_BIN = find_executable("NZVM")

siteprop = Site()
siteprop.vs30 = 500



def load_msgs(args, logger: Logger = qclogging.get_basic_logger()):
    # returns list of appropriate srf metadata
    msgs = []
    faults = set()

    # add task for every info file
    rel_files = glob(args.rel_glob)

    # add task for every info file

    for rel in rel_files:
        if not os.path.exists(rel):
            logger.info("REL csv file not found: {}".format(rel))
            continue

        # name is unique and based on basename
        name = os.path.splitext(os.path.basename(rel))[0].split("_")[0]

        if name in faults:
            continue
        logger.debug("Found first REL file for {}".format(name))
        faults.add(name)

        a=pd.read_csv(rel)

        rake= a['rake'].loc[0]


        num_subfaults = a['plane_count'].loc[0]
        corners=[]
        for i in range(num_subfaults):
            latitude = a['clat_subfault_{}'.format(i)].loc[0]
            longitude = a['clon_subfault_{}'.format(i)].loc[0]
            flen = a['length_subfault_{}'.format(i)].loc[0]
            fwid = a['width_subfault_{}'.format(i)].loc[0]
            dip = a['dip_subfault_{}'.format(i)].loc[0]
            strike = a['strike_subfault_{}'.format(i)].loc[0]
            corners.append(get_corners(latitude, longitude, flen, fwid, dip, strike))

        try:
            mag = mom2mag(sum(map(mag2mom, a["magnitude"].loc[0])))
        except TypeError:
            mag = a["magnitude"].loc[0]
        dtop = a["dtop"].loc[0]
        dbottom = a["dbottom"].loc[0]

        msgs.append(
            (
                args,
                {
                    "name": name,
                    "dip": a["dip"].loc[0],
                    "rake": rake,
                    "dbottom": dbottom,
                    "corners": np.array(corners),
                    "mag": mag,
                    "hdepth":  dtop + 0.5 * (dbottom - dtop),
                },
                qclogging.get_realisation_logger(logger, name).name,
                False,
            )
        )
    print(msgs)
    return msgs


def load_args(logger: Logger = qclogging.get_basic_logger()):
    parser = ArgumentParser()
    arg = parser.add_argument

    arg("rel_glob", help="rel csv file selection expression. eg: Srf/*.csv")


    arg(
        "-o", "--out-dir", help="directory to place outputs", default="VMs"
    )
    arg(
        "--pgv",
        help="max PGV at velocity model perimiter (estimated, cm/s)",
        type=float,
        default=5.0,
    )
    arg("--hh", help="velocity model grid spacing (km)", type=float, default=0.4)
    arg(
        "--dt",
        help="timestep to estimate simulation duration (s)",
        type=float,
        default=0.005,
    )
    arg(
        "--space-land",
        help="min space between VM edge and land (km)",
        type=float,
        default=5.0,
    )
    arg(
        "--space-srf",
        help="min space between VM edge and SRF (km)",
        type=float,
        default=15.0,
    )
    arg("--min-vs", help="for nzvm gen and flo (km/s)", type=float, default=0.5)
    arg("-n", "--nproc", help="number of processes", type=int, default=1)
    arg(
        "-t",
        "--vm_threads",
        "--threads",
        help="number of threads for the VM generation",
        type=int,
        default=1,
    )
    arg("--novm", help="only generate parameters", action="store_true")
    arg("--vm-version", help="velocity model version to generate", default="1.65")
    arg(
        "--vm-topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )
    arg("--selection", help="also generate NHM selection file", action="store_true")
    arg(
        "--no_optimise",
        help="Don't try and optimise the vm if it is off shore. Removes dependency on having GMT coastline data",
        action="store_true",
    )
    arg(
        "--min-rjb",
        help="Specify a minimum horizontal distance (in km) for the VM to span from the fault"
             " - invalid VMs will still not be generated",
        default=0,
    )
    arg(
        "--ds-multiplier",
        help="Sets the DS multiplier for setting the sim-duration. Validation runs default to 1.2. Cybershake runs"
             "should manually set it to 0.75",
        default=1.2,
        type=float,
    )
    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)



    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("rel2vm_params")
    qclogging.add_general_file_handler(
        logger, os.path.join(os.getcwd(), "rel2vm_params_log.txt")
    )
    args = load_args(logger=logger)

    args.novm=True


    logger.debug(
        "Loading Srf/*.csv"
    )
    msg_list = load_msgs(args, logger=logger)
    if len(msg_list) == 0:
        message = "Found nothing to do, exiting."
        logger.log(qclogging.NOPRINTCRITICAL, message)
        sys.exit(message)

    # prepare to run
    if not os.path.isdir(args.out_dir):
        logger.debug(
            "Output directory {} does not exist. Creating it now.".format(args.out_dir)
        )
        os.makedirs(args.out_dir)

    # distribute work
    logger.debug(
        "Number of processes to use given as {}, number of tasks: {}.".format(
            args.nproc, len(msg_list)
        )
    )
    p = Pool(processes=args.nproc)
    reports = p.starmap(create_vm, msg_list)

    # store summary
    store_summary(os.path.join(args.out_dir, "rel2vm_params_info.csv"), reports, logger=logger)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)

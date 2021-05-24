#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

most basic example (set out_dir with -o or --out-dir):
mpirun -n 3 vm_params2vm.py ./vm_params.yml Hossack

HELP:
./vm2_params2vm.py -h
"""

from argparse import ArgumentParser

from distutils.spawn import find_executable

from logging import Logger

from multiprocessing import Pool

import os
import platform
import subprocess
from shutil import rmtree, copyfile

import sys
from tempfile import mkdtemp
import yaml

from qcore import constants, qclogging

from srfinfo2vm import get_env_dirs, gen_vm, save_nzvm_cfg, store_summary, faultprop, NZVM_BIN


def load_vm_params(vm_params_file):
    with open(vm_params_file, 'r') as f:
        vm_params_dict = yaml.load(f, Loader=yaml.SafeLoader)

    return vm_params_dict


# does both vm_params and vm
def create_vm(args, name, logger_name: str = "vm_params2vm"):
    # temp directory for current process

    logger = qclogging.get_realisation_logger(
        qclogging.get_logger(logger_name), name
    )
    ptemp = mkdtemp(prefix="_tmp_%s_" % (name), dir=args.out_dir)

    vm_params_dict = args.vm_params_dict
    (out_vm_dir, vm_working_dir, nzvm_cfg, vm_params_path) = get_env_dirs(args, name, ptemp, mkdir=False)
    save_nzvm_cfg(nzvm_cfg, vm_params_dict, vm_working_dir, logger=logger)

    if vm_params_dict is not None:
        # run the actual generation
        copyfile(args.vm_params, os.path.join(ptemp, os.path.basename(args.vm_params)))
        gen_vm(args, name, vm_params_dict, faultprop.Mw, ptemp, logger=logger)
    else:
        logger.debug("At least one dimension was 0. Not generating VM")
    # plot results
    vm_params_dict["vm_dir"] = os.path.join(args.out_dir, name)

    # working dir cleanup, return info about VM
    logger.debug("Cleaning up temp directory {}".format(ptemp))
    rmtree(ptemp)
    return vm_params_dict


def load_msgs_vm_params(args, logger: Logger = qclogging.get_basic_logger()):
    msgs = []
    logger.debug("Checking VM params file: {}".format(args.vm_params))

    if not os.path.exists(args.vm_params):
        raise FileNotFoundError(args.vm_params)

    args.vm_params_dict = load_vm_params(args.vm_params)

    msgs.append(
        (
            args,
            args.name,
            qclogging.get_realisation_logger(logger, args.name).name,
        )
    )
    return msgs


def load_args(logger: Logger = qclogging.get_basic_logger()):
    parser = ArgumentParser()
    arg = parser.add_argument

    arg(
        "vm_params",
        help="path to vm_params.yaml",

    )
    arg(
        "name",
        help="Name of the model"
    )
    arg(
        "-o", "--out-dir", help="directory to place outputs", default="VMs"
    )

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

    args = parser.parse_args()
    args.out_dir = os.path.abspath(args.out_dir)

    if not args.novm:
        if NZVM_BIN is None:
            message = """NZVM binary not in PATH
You can compile it from here: https://github.com/ucgmsim/Velocity-Model
Then add it to the PATH by either:
sudo cp /location/to/Velocity-Model/NZVM /usr/bin/
or:
export PATH=$PATH:/location/to/Velocity-Model"""
            logger.log(qclogging.NOPRINTERROR, message)
            sys.exit(message)

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("vm_params2vm")
    qclogging.add_general_file_handler(
        logger, os.path.join(os.getcwd(), "vm_params2vm_log.txt")
    )
    args = load_args(logger=logger)

    # prepare to run
    if not os.path.isdir(args.out_dir):
        logger.debug(
            "Output directory {} does not exist. Creating it now.".format(args.out_dir)
        )
        os.makedirs(args.out_dir)

    msg_list = load_msgs_vm_params(args, logger=logger)

    # distribute work
    logger.debug(
        "Number of processes to use given as {}, number of tasks: {}.".format(
            args.nproc, len(msg_list)
        )
    )
    p = Pool(processes=args.nproc)
    reports = p.starmap(create_vm, msg_list)

    # store summary
    store_summary(os.path.join(args.out_dir, "vminfo.csv"), reports, logger=logger)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)

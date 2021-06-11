#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

Most basic example:

mpirun -n 3 vm_params2vm.py Hossack ./vm_params.yml ~/Data/VMs
(The root directory of VMs should be given. ~/Data/VMs/Hossack will be auto-created)

HELP:
./vm2_params2vm.py -h
"""

from argparse import ArgumentParser
from distutils.spawn import find_executable

from logging import Logger
from multiprocessing import Pool
import os
import platform

from shutil import rmtree, copyfile, move
import subprocess
from subprocess import Popen
import sys

from tempfile import mkdtemp
import yaml

from qcore import constants, qclogging
from qcore.validate_vm import validate_vm

from common import store_summary, temp_paths
from gen_coords import gen_coords

NZVM_BIN = find_executable("NZVM")
NZVM_BIN = "/home/jam335/code/Velocity-Model/NZVM"


def gen_vm(
    out_dir,
    ptemp,
    vm_threads=1,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Generates VM and relocate the output files from ptemp to out_dir
    Also generates coordinates files and validates.

    Parameters
    ----------
    out_dir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    ptemp : Temporary work directory (eg. .../Data/VMs/Hossack/_tmp_Hossack_nho8ui41)
    vm_threads : Number of threads to run NZVM binary with if OpenMP-enabled
    logger :
    """
    os.makedirs(out_dir, exist_ok=True)
    vm_working_dir, nzvm_cfg, vm_params_path = temp_paths(ptemp)

    # NZVM won't find resources if WD is not NZVM dir, stdout not MPROC friendly
    with open(os.path.join(ptemp, "NZVM.out"), "w") as logfile:
        logger.debug("Running NZVM binary")
        nzvm_env = os.environ.copy()
        nzvm_env["OMP_NUM_THREADS"] = str(vm_threads)
        nzvm_exe = Popen(
            [NZVM_BIN, nzvm_cfg],
            cwd=os.path.dirname(NZVM_BIN),
            stdout=logfile,
            env=nzvm_env,
        )
        nzvm_exe.communicate()
    logger.debug("Moving VM files to vm directory")
    # fix up directory contents
    logger.debug("Removing Log and Velocity_Model directories")
    logger.debug("Moving nzvm config and vm_params yaml to vm directory")
    for f in [
        os.path.join(vm_working_dir, "Velocity_Model", vm3dfile)
        for vm3dfile in ["rho3dfile.d", "vp3dfile.p", "vs3dfile.s"]
    ]:
        move(f, os.path.join(out_dir, os.path.basename(f)))  # may overwrite

    logger.debug("Removing Log and Velocity_Model directories")
    for f in [
        os.path.join(vm_working_dir, "Log", "VeloModCorners.txt"),
        nzvm_cfg,
        vm_params_path,
    ]:
        move(f, os.path.join(out_dir, os.path.basename(f)))  # may overwrite

    rmtree(vm_working_dir)
    # create model_coords, model_bounds etc...
    logger.debug("Generating coords")
    gen_coords(vm_dir=out_dir)
    # validate
    logger.debug("Validating vm")
    success, message = validate_vm(out_dir)
    if success:
        vm_check_str = f"VM check passed: {out_dir}"
        logger.debug(vm_check_str)
        print(vm_check_str, file=sys.stderr)
    else:
        logger.log(
            qclogging.NOPRINTCRITICAL,
            "VM check for {} failed: {}".format(out_dir, message),
        )
        print(f"VM check BAD: {message}", file=sys.stderr)


def save_nzvm_cfg(
    nzvm_cfg_path,
    vm_params_dict,
    vm_dir,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Saves nzvm.cfg from vm_params_dict.
    NZVM binary needs to be run from its installed directory with a path to nzvm.cfg
    Its outputs are to be placed in OUTPUT_DIR, which should not pre-exist. OUTPUT_DIR should be in absolute path

    Parameters
    ----------
    nzvm_cfg_path : Path to nzvm.cfg (initially placed in ptemp, and relocated after NZVM run is successful)
    vm_params_dict : Dictionary containing relevant data
    vm_dir : Directory where NZVM binary outputs are placed (cannot exist)
    logger :
    """
    if nzvm_cfg_path is not None:
        with open(nzvm_cfg_path, "w") as vmd:
            vmd.write(
                "\n".join(
                    [
                        "CALL_TYPE=GENERATE_VELOCITY_MOD",
                        "MODEL_VERSION=%s" % (vm_params_dict["model_version"]),
                        "OUTPUT_DIR=%s" % (vm_dir),
                        "ORIGIN_LAT=%s" % (vm_params_dict["MODEL_LAT"]),
                        "ORIGIN_LON=%s" % (vm_params_dict["MODEL_LON"]),
                        "ORIGIN_ROT=%s" % (vm_params_dict["MODEL_ROT"]),
                        "EXTENT_X=%s" % (vm_params_dict["extent_x"]),
                        "EXTENT_Y=%s" % (vm_params_dict["extent_y"]),
                        "EXTENT_ZMAX=%s" % (vm_params_dict["extent_zmax"]),
                        "EXTENT_ZMIN=%s" % (vm_params_dict["extent_zmin"]),
                        "EXTENT_Z_SPACING=%s" % (vm_params_dict["hh"]),
                        "EXTENT_LATLON_SPACING=%s" % (vm_params_dict["hh"]),
                        "MIN_VS=%s" % (vm_params_dict["min_vs"]),
                        "TOPO_TYPE=%s\n" % (vm_params_dict["topo_type"]),
                    ]
                )
            )
        logger.debug("Saved nzvm config file to {}".format(nzvm_cfg_path))
    else:
        logger.debug("Path to nzvm.conf should be specified")


def main(
    name, vm_params_path, out_dir, vm_threads=1, logger_name: str = "vm_params2vm"
):
    """
    Orchestrates conversion from vm_params.yaml to VM output

    Parameters
    ----------
    name : name of the fault
    vm_params_path : Path to the vm_params.yaml
    out_dir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    vm_threads : Number of threads to run NZVM binary with if OpenMP-enabled
    logger_name : Name of the logger

    Returns
    -------
    vm_params_dict : Used for store_summary
    """
    logger = qclogging.get_realisation_logger(qclogging.get_logger(logger_name), name)

    # temp directory for current process
    ptemp = mkdtemp(prefix="_tmp_%s_" % (name), dir=out_dir)

    qclogging.add_general_file_handler(
        logger, os.path.join(ptemp, "vm_params2vm_{}_log.txt".format(name))
    )

    (vm_working_dir, nzvm_cfg_path, _) = temp_paths(ptemp)

    with open(vm_params_path, "r") as f:
        vm_params_dict = yaml.load(f, Loader=yaml.SafeLoader)

    if vm_params_dict is not None:
        # saves nzvm.cfg
        save_nzvm_cfg(nzvm_cfg_path, vm_params_dict, vm_working_dir, logger=logger)
        # makes a copy of vm_params.yaml
        copyfile(vm_params_path, os.path.join(ptemp, os.path.basename(vm_params_path)))
        # run the actual generation
        gen_vm(out_dir, ptemp, vm_threads=vm_threads, logger=logger)
    else:
        logger.debug("vm_params.yaml has no data for VM generation")

    vm_params_dict["out_dir"] = out_dir

    # working dir cleanup, return info about VM
    logger.debug("Cleaning up temp directory {}".format(ptemp))
    rmtree(ptemp)
    return vm_params_dict


def load_msgs_vm_params(args, logger: Logger = qclogging.get_basic_logger()):

    msgs = []
    logger.debug("Checking VM params file: {}".format(args.vm_params_path))

    if not os.path.exists(args.vm_params_path):
        raise FileNotFoundError(args.vm_params_path)

    msgs.append(
        (
            args.name,
            args.vm_params_path,
            args.out_dir,
            args.vm_threads,
            qclogging.get_realisation_logger(logger, args.name).name,
        )
    )
    return msgs


def load_args(logger: Logger = qclogging.get_basic_logger()):

    parser = ArgumentParser()
    arg = parser.add_argument

    arg("name", help="Name of the fault")

    arg(
        "vm_params_path",
        help="path to vm_params.yaml",
    )

    arg("out_dir", help="root directory to place VM files eg. Data/VMs", type=os.path.abspath)

    arg("-n", "--nproc", help="number of processes", type=int, default=1)
    arg(
        "-t",
        "--vm_threads",
        "--threads",
        help="number of threads for the VM generation",
        type=int,
        default=1,
    )
    arg("--vm-version", help="velocity model version to generate", default="1.65")
    arg(
        "--vm-topo",
        help="topo_type parameter for velocity model generation",
        default="BULLDOZED",
    )

    args = parser.parse_args()
    args.vm_params_path = os.path.abspath(args.vm_params_path)
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
            "VM output directory {} does not exist. Creating it now.".format(
                args.out_dir
            )
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
    reports = p.starmap(main, msg_list)

    # store summary
    store_summary(
        os.path.join(args.out_dir, "vm_params2vm_info.csv"), reports, logger=logger
    )

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)

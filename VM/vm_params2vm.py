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
import csv
from distutils.spawn import find_executable

from logging import Logger
from multiprocessing import Pool
import os
import platform

from shutil import rmtree, copyfile, move
import subprocess
from subprocess import Popen
import sys


from tempfile import TemporaryDirectory
import yaml

from qcore import constants, qclogging
from qcore.validate_vm import validate_vm

from gen_coords import gen_coords

NZVM_BIN = find_executable("NZVM")

def temp_paths(ptemp):
    vm_working_dir = os.path.join(ptemp, "output")
    nzvm_cfg = os.path.join(ptemp, "nzvm.cfg")
    vm_params_path = os.path.join(ptemp, "vm_params.yaml")
    return (vm_working_dir, nzvm_cfg, vm_params_path)


def store_summary(table, info_store, logger: Logger = qclogging.get_basic_logger()):
    # initialise table file
    logger.debug("Saving summary")
    keys = info_store[0].keys()
    with open(table, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        writer.writerow(keys)
        for i in info_store:
            writer.writerow([i[key] for key in keys])
    logger.info("Saved summary to {}".format(table))




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

    logger.debug("Moving output files to vm directory")
    files_to_move = [
        os.path.join(vm_working_dir, "Velocity_Model", vm3dfile)
        for vm3dfile in ["rho3dfile.d", "vp3dfile.p", "vs3dfile.s"]
    ] +  [
        os.path.join(vm_working_dir, "Log", "VeloModCorners.txt"),
        nzvm_cfg,
        vm_params_path,
    ]

    for f in files_to_move:
        move(f, os.path.join(out_dir, os.path.basename(f)))  # may overwrite
    logger.debug("Removing Log and Velocity_Model directories")

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
                        f"MODEL_VERSION={vm_params_dict['model_version']}",
                        f"OUTPUT_DIR={vm_dir}",
                        f"ORIGIN_LAT={vm_params_dict['MODEL_LAT']}",
                        f"ORIGIN_LON={vm_params_dict['MODEL_LON']}",
                        f"ORIGIN_ROT={vm_params_dict['MODEL_ROT']}",
                        f"EXTENT_X={vm_params_dict['extent_x']}",
                        f"EXTENT_Y={vm_params_dict['extent_y']}",
                        f"EXTENT_ZMAX={vm_params_dict['extent_zmax']}",
                        f"EXTENT_ZMIN={vm_params_dict['extent_zmin']}",
                        f"EXTENT_Z_SPACING={vm_params_dict['hh']}",
                        f"EXTENT_LATLON_SPACING={vm_params_dict['hh']}",
                        f"MIN_VS={vm_params_dict['min_vs']}",
                        f"TOPO_TYPE={vm_params_dict['topo_type']}\n",
                    ]
                )
            )
        logger.debug("Saved nzvm config file to {}".format(nzvm_cfg_path))
    else:
        logger.debug("Path to nzvm.conf should be specified")


def main(
    name, vm_params_path, out_dir, vm_threads=1, logger: Logger = qclogging.get_basic_logger()
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


    # temp directory for current process
    with TemporaryDirectory(prefix=f"_tmp_{name}_", dir=out_dir) as ptemp:

        qclogging.add_general_file_handler(
            logger, os.path.join(out_dir, "vm_params2vm_{}_log.txt".format(name))
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


    return vm_params_dict



def load_args(logger: Logger = qclogging.get_basic_logger()):

    parser = ArgumentParser()
    arg = parser.add_argument

    arg("name", help="Name of the fault")

    arg(
        "vm_params_path",
        help="path to vm_params.yaml",
    )

    arg(
        "-o",
        "--out_dir",
        help="output directory to place VM files "
             "(if not specified, the same path as rel_file is used",
        default=None,
    )


    arg(
        "-t",
        "--vm_threads",
        "--threads",
        help="number of threads for the VM generation",
        type=int,
        default=1,
    )


    args = parser.parse_args()
    args.vm_params_path = os.path.abspath(args.vm_params_path)

    if args.out_dir is None:
        args.out_dir = os.path.dirname(os.path.dirname(args.vm_params_path))
    args.out_dir = os.path.abspath(args.out_dir)

    return args


if __name__ == "__main__":

    logger = qclogging.get_logger("vm_params2vm")
    qclogging.add_general_file_handler(
        logger, os.path.join(os.getcwd(), "vm_params2vm_log.txt")
    )
    args = load_args(logger=logger)

    # prepare to run
    logger.debug("Checking VM params file: {}".format(args.vm_params_path))
    if not os.path.exists(args.vm_params_path):
        raise FileNotFoundError(args.vm_params_path)

    logger.debug("Checking/Creating output directory :{}".format(args.out_dir))
    os.makedirs(args.out_dir, exist_ok=True)

    main(args.name, args.vm_params_path, args.out_dir,args.vm_threads,logger=logger)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)

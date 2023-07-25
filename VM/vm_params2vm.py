#!/usr/bin/env python3
"""
REQUIREMENTS:
PATH contains NZVM binary (github:ucgmsim/Velocity-Model/NZVM)

Most basic example:

python vm_params2vm.py Hossack ./vm_params.yaml -o ~/Data/VMs/Hossack

HELP:
python vm2_params2vm.py -h
"""

from argparse import ArgumentParser
from distutils.spawn import find_executable

from logging import Logger
import os
from pathlib import Path
from shutil import rmtree, copyfile, move
from subprocess import Popen
import sys


from tempfile import TemporaryDirectory,mkdtemp
import yaml

from qcore import qclogging
from qcore.validate_vm import validate_vm_files

from gen_coords import gen_coords

NZVM_BIN = find_executable("NZVM")


def gen_vm(
    outdir: Path,
    temp_dir: Path,
    vm_working_dir: Path,
    nzvm_cfg_path: Path,
    vm_params_path: Path,
    vm_threads=1,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Generates VM and relocate the output files from temp_dir to outdir
    Also generates coordinates files and validates.

    Parameters
    ----------
    outdir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    temp_dir : Temporary work directory (eg. .../Data/VMs/Hossack/_tmp_Hossack_nho8ui41)
    vm_working_dir: working directory that NZVM will be generating its output/temporary files
    nzvm_cfg_path: the path to nzvm.cfg
    vm_params_path: the path to vm_params.yaml
    vm_threads : Number of threads to run NZVM binary with if OpenMP-enabled
    logger :
    """
    os.makedirs(outdir, exist_ok=True)
    # NZVM won't find resources if WD is not NZVM dir, stdout not MPROC friendly
    with open(temp_dir / "NZVM.out", "w") as logfile:
        logger.debug("Running NZVM binary")
        nzvm_env = os.environ.copy()
        nzvm_env["OMP_NUM_THREADS"] = str(vm_threads)
        nzvm_exe = Popen(
            [NZVM_BIN, nzvm_cfg_path],
            cwd=Path(NZVM_BIN).parent,
            stdout=logfile,
            env=nzvm_env,
        )
        nzvm_exe.communicate()
    print(f"{NZVM_BIN} {nzvm_cfg_path}")

    with open(temp_dir /"NZVM.out", "r") as logfile:
        lines = logfile.readlines()
        logger.debug("\n".join(lines))
    logger.debug([x.name for x in temp_dir.glob("*")])

    logger.debug("Moving output files to vm directory")
    files_to_move = [
        vm_working_dir / "Velocity_Model" / vm3dfile
        for vm3dfile in ["rho3dfile.d", "vp3dfile.p", "vs3dfile.s", "in_basin_mask.b"]
    ] + [nzvm_cfg_path] + list((vm_working_dir / "Log").glob("VeloModCorners*.txt")) # due to NZVM's MPI support, it is now VeloModCorners-%d

    for f in files_to_move:
        move(f, outdir / f.name)  # may overwrite

    vm_params_path_dest = outdir / vm_params_path.name
    if vm_params_path != vm_params_path_dest:
        if vm_params_path_dest.is_file():
            vm_params_path_dest.unlink()
        copyfile(vm_params_path, vm_params_path_dest)
    # otherwise, vm_params.yaml is already in outdir

    logger.debug("Removing Log and Velocity_Model directories")

    rmtree(vm_working_dir)
    # create model_coords, model_bounds etc...
    logger.debug("Generating coords")
    gen_coords(vm_dir=outdir)
    # validate
    logger.debug("Validating vm")
    success, message = validate_vm_files(outdir)
    if success:
        vm_check_str = f"VM check passed: {outdir}"
        logger.debug(vm_check_str)
        print(vm_check_str, file=sys.stderr)
    else:
        logger.log(
            qclogging.NOPRINTCRITICAL, f"VM check for {outdir} failed: {message}"
        )
        print(f"VM check BAD: {message}", file=sys.stderr)


def save_nzvm_cfg(
    nzvm_cfg_path: Path,
    vm_params_dict: dict,
    vm_dir: Path,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Saves nzvm.cfg from vm_params_dict.
    NZVM binary needs to be run from its installed directory with a path to nzvm.cfg
    Its outputs are to be placed in OUTPUT_DIR, which should not pre-exist. OUTPUT_DIR should be in absolute path

    Parameters
    ----------
    nzvm_cfg_path : Path to nzvm.cfg (initially placed in temp_dir, and relocated after NZVM run is successful)
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
        logger.debug(f"Saved nzvm config file to {nzvm_cfg_path}")
    else:
        logger.debug("Path to nzvm.conf should be specified")


def main(
    name: str,
    vm_params_path: Path,
    outdir: Path,
    vm_threads=1,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Orchestrates conversion from vm_params.yaml to VM output

    Parameters
    ----------
    name : name of the fault
    vm_params_path : Path to the vm_params.yaml
    outdir : Directory where output files are eventually saved (eg. ..../Data/VMs/Hossack)
    vm_threads : Number of threads to run NZVM binary with if OpenMP-enabled
    logger :

    """

    # temp directory for current process
#    with TemporaryDirectory(prefix=f"_tmp_{name}_", dir=outdir) as temp_dir:
#    tempd = TemporaryDirectory(prefix=f"_tmp_{name}_", dir=outdir)
    tempd = mkdtemp(prefix=f"_tmp_{name}_",dir=outdir)
    if True:
        temp_dir = Path(tempd)
        temp_dir.exists()
        qclogging.add_general_file_handler(
            logger, outdir / f"vm_params2vm_{name}_log.txt"
        )
        logger.debug(f"{temp_dir} created")

        vm_working_dir = temp_dir / "output"
        nzvm_cfg_path = temp_dir / "nzvm.cfg"

        with open(vm_params_path, "r") as f:
            vm_params_dict = yaml.load(f, Loader=yaml.SafeLoader)

        if vm_params_dict is not None:
            # saves nzvm.cfg
            save_nzvm_cfg(nzvm_cfg_path, vm_params_dict, vm_working_dir, logger=logger)

            # run the actual generation
            gen_vm(
                outdir,
                temp_dir,
                vm_working_dir,
                nzvm_cfg_path,
                vm_params_path,
                vm_threads=vm_threads,
                logger=logger,
            )
        else:
            logger.debug("vm_params.yaml has no data for VM generation")
#    tempd.cleanup()

def load_args(logger: Logger = qclogging.get_basic_logger()):
    """
    Unpacks arguments and does basic checks

    Parameters
    ----------
    logger :

    Returns
    -------
    Processed arguments
    """

    parser = ArgumentParser()
    arg = parser.add_argument

    arg("name", help="Name of the fault")

    arg("vm_params_path", help="path to vm_params.yaml")

    arg(
        "-o",
        "--outdir",
        help="output directory to place VM files "
        "(if not specified, the directory containing vm_params.yaml is used",
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
    args.vm_params_path = Path(args.vm_params_path).resolve()

    logger.debug(f"Checking VM params file: {args.vm_params_path}")
    if not args.vm_params_path.exists():
        parser.error(f"File is not present: {args.vm_params_path}")

    if args.outdir is None:
        args.outdir = args.vm_params_path.parent
    args.outdir = Path(args.outdir).resolve()

    if (
        args.outdir / args.vm_params_path.name
    ).exists() and args.outdir != args.vm_params_path.parent:
        parser.error(
            f"{args.vm_params_path.name} already exists in {args.outdir}. Aborting."
        )

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("vm_params2vm")
    qclogging.add_general_file_handler(logger, Path.cwd() / "vm_params2vm_log.txt")
    args = load_args(logger=logger)

    # prepare to run

    logger.debug(f"Checking/Creating output directory :{args.outdir}")
    os.makedirs(args.outdir, exist_ok=True)

    main(args.name, args.vm_params_path, args.outdir, args.vm_threads, logger=logger)

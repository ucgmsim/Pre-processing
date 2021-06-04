"""Generate Velocity Models for all events in a fault selection file
This script generates velocity models for all events in a fault selection file.
It requires that the first/only realisation for each event/fault has had a .info file generated.
All  command line arguments are the same as srfinfo2vm
"""
import argparse
from h5py import File as h5open
from multiprocessing import Pool
import os
from os.path import abspath, join, exists
import platform
import subprocess
import sys

from qcore import simulation_structure, qclogging
from qcore.formats import load_fault_selection_file as load_fsf

from createSRF import mag2mom, mom2mag
from VM.deprecated.srfinfo2vm import create_vm, NZVM_BIN, store_summary


NO_NZVM_MESSAGE = """NZVM binary not in PATH
You can compile it from here: https://github.com/ucgmsim/Velocity-Model
Then add it to the PATH by either:
sudo cp /location/to/Velocity-Model/NZVM /usr/bin/
or:
export PATH=$PATH:/location/to/Velocity-Model"""


def load_args(logger=qclogging.get_basic_logger()):
    """Loads all the arguments required and verifys them before returning to execution"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cs_root", type=abspath, help="Path to the root of the simulation directory"
    )
    parser.add_argument(
        "fsf_path", type=abspath, help="Path to the fault selection file"
    )
    parser.add_argument(
        "checkpointing",
        type="store_true",
        help="Prevents regeneration of existing VMs. Checks for the presence of the vm_params.yaml file.",
    )
    parser.add_argument(
        "--pgv",
        help="Max PGV at velocity model perimiter (estimated, cm/s)",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--hh",
        help="Velocity model grid spacing (km). Default value: 0.4",
        type=float,
        default=0.4,
    )
    parser.add_argument(
        "--dt",
        help="Timestep to estimate simulation duration (s). Default value: 0.02, compatible with all current gmsim versions (down to 100m)",
        type=float,
        default=0.02,
    )
    parser.add_argument(
        "--space-land",
        help="Minimum space between VM edge and land (km). Default value: 5",
        type=float,
        default=5.0,
    )
    parser.add_argument(
        "--space-srf",
        help="Minimum space between VM edge and SRF (km). Default value: 15",
        type=float,
        default=15.0,
    )
    parser.add_argument(
        "--min-vs",
        help="for nzvm gen and flo (km/s). Default value: 0.5",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "-n",
        "--nproc",
        help="Number of processes. Default value: 1",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-t",
        "--vm_threads",
        "--threads",
        help="Number of threads for each NVZM instance to use. Default value: 1",
        type=int,
        default=1,
    )
    parser.add_argument("--novm", help="Only generate parameters", action="store_true")
    parser.add_argument(
        "--vm-version",
        help="Velocity model version to generate. Default value: 1.65",
        default="1.65",
    )
    parser.add_argument(
        "--vm-topo",
        help="Topo_type parameter for velocity model generation. Default value: BULLDOZED",
        default="BULLDOZED",
    )
    parser.add_argument(
        "--no_optimise",
        help="Don't try and optimise the vm if it is off shore. Removes dependency on having GMT coastline data",
        action="store_true",
    )
    parser.add_argument(
        "--min-rjb",
        help="Specify a minimum horizontal distance (in km) for the VM to span from the fault - invalid VMs will still not be generated. Default value: 0",
        default=0,
    )
    parser.add_argument(
        "--ds-multiplier",
        help="Sets the DS multiplier for setting the sim-duration. Validation runs default to 1.2. Cybershake runs"
        "should manually set it to 0.75",
        default=1.2,
        type=float,
    )
    args = parser.parse_args()
    if not args.novm and NZVM_BIN is None:
        sys.exit(NO_NZVM_MESSAGE)
    return args


def generate_tasks(args, logger=qclogging.get_basic_logger()):
    """
    Generates the tasks to run
    Returns a list of tasks that can be sent to create_vm with starmap
    :param args: The args object returned by ArgumentParser
    """
    tasks = []
    faults = load_fsf(args.fsf_path)
    checkpointing = args.checkpointing
    for f, c in faults.items():
        if checkpointing and exists(
            join(
                simulation_structure.get_fault_VM_dir(args.cs_root, f), "vm_params.yaml"
            )
        ):
            logger.info(f"VM params file found for event/fault {f}, not generating")
            continue

        # Grab either the only realisation, or the first
        # Change this to the unperturbated info when properly implemented
        if c == 1:
            rel = f
        else:
            rel = simulation_structure.get_realisation_name(f, 1)
        info_file = join(
            simulation_structure.get_sources_dir(args.cs_root),
            simulation_structure.get_srf_info_location(rel),
        )
        with h5open(info_file, "r") as h:
            a = h.attrs
            try:
                rake = a["rake"][0][0]
            except (TypeError, IndexError):
                rake = a["rake"]
            try:
                rake = a["rake"][0][0]
            except (TypeError, IndexError):
                rake = a["rake"]
            try:
                mag = mom2mag(sum(map(mag2mom, a["mag"])))
            except TypeError:
                mag = a["mag"]
            tasks.append(
                (
                    args,
                    {
                        "name": f,
                        "dip": a["dip"][0],
                        "rake": rake,
                        "dbottom": a["dbottom"][0],
                        "corners": a["corners"],
                        "mag": mag,
                        "hdepth": a["hdepth"],
                    },
                    qclogging.get_realisation_logger(logger, f).name,
                )
            )
    return tasks


def main():
    """Main function when running this as a script. Prevents globals"""
    logger = qclogging.get_logger("srfinfo2vm")
    args = load_args(logger=logger)

    qclogging.add_general_file_handler(logger, join(args.cs_root, "srfinfo2vm_log.txt"))

    out_dir = simulation_structure.get_VM_dir(args.cs_root)
    args.out_dir = out_dir

    task_list = generate_tasks(args, logger=logger)
    if len(task_list) == 0:
        logger.critical("Found nothing to do, exiting")
        return

    os.makedirs(out_dir, exist_ok=True)

    p = Pool(processes=args.nproc)
    reports = p.starmap(create_vm, task_list)

    store_summary(join(out_dir, "vminfo.csv"), reports, logger=logger)

    # Hack to fix VM generation permission issue
    hostname = platform.node()
    if hostname.startswith(("maui", "mahuika", "wb", "ni")):  # Checks if is on the HPCF
        logger.debug("HPC detected, changing folder ownership and permissions")
        permission_cmd = ["chmod", "g+rwXs", "-R", args.out_dir]
        subprocess.call(permission_cmd)
        group_cmd = ["chgrp", constants.DEFAULT_ACCOUNT, "-R", args.out_dir]
        subprocess.call(group_cmd)


if __name__ == "__main__":
    main()

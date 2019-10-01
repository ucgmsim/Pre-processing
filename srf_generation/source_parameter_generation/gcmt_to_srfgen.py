#! /usr/bin/python3

import argparse
from logging import Logger
from multiprocessing import pool
from os import makedirs
from os.path import abspath, isfile, dirname, join
from typing import Callable

import pandas as pd
import yaml

from qcore.formats import load_fault_selection_file
from qcore.simulation_structure import get_realisation_name, get_srf_path
from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_realisation_logger,
)

# from SrfGen.source_parameter_generation.uncertainties.versions import PERTURBATORS
from srf_generation.source_parameter_generation.uncertainties.common import (
    GCMT_PARAM_NAMES,
    GCMT_PARAMS,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
)


def load_args(primary_logger: Logger):
    parser = argparse.ArgumentParser()

    parser.add_argument("version", type=str)
    parser.add_argument("fault_selection_file")
    parser.add_argument("gcmt_file")
    parser.add_argument("-n", "--n_processes", default=1)
    parser.add_argument("-c", "--cybershake_root", type=str, default=".")

    args = parser.parse_args()

    args.fault_selection_file = abspath(args.fault_selection_file)
    args.gcmt_file = abspath(args.gcmt_file)

    errors = []

    if not isfile(args.station_file):
        errors.append(f"Specified station file not found: {args.station_file}")
    if not isfile(args.gcmt_file):
        errors.append(f"Specified gcmt file not found: {args.gcmt_file}")

    if errors:
        message = (
            "At least one error was detected when verifying arguments:\n"
            + "\n".join(errors)
        )
        raise ValueError(message)

    makedirs(args.cybershake_root, exist_ok=True)

    log_file = join(args.cybershake_root, primary_logger.name + "_log.txt")
    add_general_file_handler(primary_logger, log_file)

    primary_logger.debug(f"Arguments parsed: {args}")

    return args


def generate_uncertainties(
    data: GCMT_PARAMS,
    realisation_count: int,
    cybershake_root: str,
    perturbation_function: Callable,
    primary_logger: Logger,
):
    fault_logger = get_realisation_logger(primary_logger, data.pid)
    fault_logger.debug(f"Fault {data.pid} had data {data}")
    for i in range(realisation_count):
        fault_logger.debug(
            f"Generating realisation {i} of {realisation_count} for fault {data.pid}"
        )
        rel_name = get_realisation_name(data.pid, i)

        fault_logger.debug("Calling perturbation function.")
        perturbed_realisation = perturbation_function(data)
        fault_logger.debug(
            f"Got results from perturbation_function: {perturbed_realisation}"
        )

        file_name = get_srf_path(cybershake_root, rel_name).replace(".srf", ".yaml")
        makedirs(dirname(file_name), exist_ok=True)
        fault_logger.debug(
            f"Created Srf directory and attempting to save perturbated source generation parameters there: {file_name}"
        )
        with open(file_name) as params_file:
            yaml.dump(perturbed_realisation, params_file)
        fault_logger.debug(
            f"Parameters saved succesfully. Continuing to next realisation if one exists."
        )


def main():

    primary_logger = get_logger("GCMT_2_srfgen_params")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    faults = load_fault_selection_file(args.fault_selection_file)
    gcmt_data = pd.read_csv(
        args.gcmt_file, usecols=(0, 2, 3, 13, 11, 4, 5, 6), names=GCMT_PARAM_NAMES
    )
    primary_logger.debug(
        f"{len(faults.keys())} faults were selected and gcmt data for {len(gcmt_data)} faults were loaded. "
        f"Running cross comparison now."
    )

    gcmt_lines = gcmt_data[gcmt_data["PublicID"].isin(faults.keys())].itertuples(
        name="GCMT_Source", index=False
    )
    missing_faults = [
        fault for fault in faults.keys() if fault not in gcmt_lines["PublicID"]
    ]
    if missing_faults:
        message = f"Some fault(s) in the fault selection were not found in the gcmt file: {missing_faults}"
        primary_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)
    primary_logger.debug(
        f"All {len(gcmt_lines)} faults specified in the fault selection file were found"
    )

    messages = [
        (
            gcmt_lines[i],
            faults[gcmt_lines[i].PublicID],
            args.cybershake_root,
            perturbation_function,
            primary_logger,
        )
        for i in range(len(faults))
    ]

    if len(messages) < args.n_processes:
        n_processes = len(messages)
        primary_logger.info(
            f"The number of processes requested to create realisations ({args.n_processes}) is greater than "
            f"the number of faults ({n_processes}), this value will be capped to the number of faults"
        )
    else:
        n_processes = args.n_processes

    primary_logger.debug(
        f"{len(messages)} messages were created to create a total of {sum([m[1] for m in messages])} realisations"
    )

    worker_pool = pool.Pool(processes=n_processes)
    worker_pool.starmap(generate_uncertainties, messages)


if __name__ == "__main__":
    main()

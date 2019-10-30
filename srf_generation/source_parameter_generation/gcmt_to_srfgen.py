#! /usr/bin/python3

import argparse
from logging import Logger
from multiprocessing import pool
from os import makedirs
from os.path import abspath, isfile, dirname, join
from typing import Callable, Union

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
    GCMT_Source,
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
    parser.add_argument("-a", "--aggregate_file", type=str)

    args = parser.parse_args()

    args.fault_selection_file = abspath(args.fault_selection_file)
    args.gcmt_file = abspath(args.gcmt_file)
    args.cybershake_root = abspath(args.cybershake_root)
    if args.aggregate_file is not None:
        args.aggregate_file = abspath(args.aggregate_file)

    errors = []

    if not isfile(args.fault_selection_file):
        errors.append(f"Specified station file not found: {args.fault_selection_file}")
    if not isfile(args.gcmt_file):
        errors.append(f"Specified gcmt file not found: {args.gcmt_file}")
    if args.aggregate_file is not None and isfile(args.aggregate_file):
        errors.append(
            f"Specified aggregation file {args.aggregate_file} already exists, please choose another file"
        )

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
    data: GCMT_Source,
    realisation_count: int,
    cybershake_root: str,
    perturbation_function: Callable,
    aggregate_file: Union[str, None],
    primary_logger_name: str,
):
    primary_logger = get_logger(name=primary_logger_name)
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

        file_name = get_srf_path(cybershake_root, rel_name).replace(".srf", ".csv")
        makedirs(dirname(file_name), exist_ok=True)
        fault_logger.debug(
            f"Created Srf directory and attempting to save perturbated source generation parameters there: {file_name}"
        )
        rel_df = pd.DataFrame(perturbed_realisation, index=[i])

        with open(file_name, "w") as params_file:
            rel_df.to_csv(file_name)

        if aggregate_file is not None:
            if not isfile(aggregate_file):
                rel_df.to_csv(aggregate_file)
            else:
                rel_df.to_csv(aggregate_file, mode="a", header=False)

        fault_logger.debug(
            f"Parameters saved succesfully. Continuing to next realisation if one exists."
        )


def main():

    primary_logger = get_logger("GCMT_2_srfgen_params")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    faults = load_fault_selection_file(args.fault_selection_file)
    types = [
        "object",
        "float64",
        "float64",
        "float64",
        "float64",
        "float64",
        "float64",
        "float64",
    ]
    dtype = {GCMT_PARAM_NAMES[i]: types[i] for i in range(len(types))}

    gcmt_data = pd.read_csv(
        args.gcmt_file,
        usecols=[
            "PublicID",
            "Latitude",
            "Longitude",
            "strike1",
            "dip1",
            "rake1",
            "Mw",
            "CD",
        ],
        dtype=dtype,
    )[["PublicID", "Latitude", "Longitude", "CD", "Mw", "strike1", "dip1", "rake1"]]
    gcmt_data.columns = GCMT_PARAM_NAMES

    primary_logger.debug(
        f"{len(faults.keys())} faults were selected and gcmt data for {len(gcmt_data)} faults were loaded. "
        f"Running cross comparison now."
    )

    gcmt_lines = list(
        gcmt_data[gcmt_data["pid"].isin(faults.keys())].itertuples(
            name=None, index=False
        )
    )

    pids = [gcmt[0] for gcmt in gcmt_lines]
    missing_faults = [fault for fault in faults.keys() if fault not in pids]

    if missing_faults:
        message = f"Some fault(s) in the fault selection were not found in the gcmt file: {missing_faults}"
        primary_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)

    primary_logger.debug(
        f"All {len(gcmt_lines)} faults specified in the fault selection file were found"
    )

    messages = [
        (
            GCMT_Source(*gcmt_lines[i]),
            faults[gcmt_lines[i][0]],
            args.cybershake_root,
            perturbation_function,
            args.aggregate_file,
            primary_logger.name,
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
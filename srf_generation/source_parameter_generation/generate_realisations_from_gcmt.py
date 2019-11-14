#! /usr/bin/python3

import argparse
from logging import Logger
from multiprocessing import pool
from os import makedirs
from os.path import abspath, isfile, join
from typing import Callable, Union, Dict, Any

import pandas as pd

from qcore.formats import load_fault_selection_file
from qcore.simulation_structure import (
    get_realisation_name,
    get_srf_path,
    get_sources_dir,
    get_realisation_VM_dir,
)
from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_realisation_logger,
)
from srf_generation.source_parameter_generation.gcmt_to_realisation import (
    GCMT_DTYPE,
    load_1d_velocity_mod,
    GCMT_FILE_COLUMNS,
    generate_realisation,
    UNPERTURBATED,
)

from srf_generation.source_parameter_generation.uncertainties.common import (
    GCMT_PARAM_NAMES,
    GCMT_Source,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
)


def load_args(primary_logger: Logger):
    parser = argparse.ArgumentParser(
        "Generates realisation files from gcmt and other input data. "
        "Intended to be used with the standard qcore simultaion structure."
    )

    parser.add_argument("fault_selection_file", type=abspath)
    parser.add_argument("gcmt_file", type=abspath)
    parser.add_argument("--version", type=str, default="unperturbated")
    parser.add_argument("--vel_mod_1d", type=abspath)
    parser.add_argument("--cybershake_root", type=abspath, default=abspath("."))
    parser.add_argument(
        "-n",
        "--n_processes",
        default=1,
        type=int,
        help="Number of processes to generate realisations with. "
        "Setting this higher than the number of faults will not have any additional effect.",
    )
    parser.add_argument("-c", "--output_directory", type=abspath, default=abspath("."))
    parser.add_argument(
        "-a",
        "--aggregate_file",
        type=abspath,
        help="A filepath to the location an aggregate file should be stored. "
        "There should not be a file already present.",
    )
    parser.add_argument(
        "--source_parameter",
        type=str,
        nargs=2,
        action="append",
        metavar=("name", "filepath"),
        help="Values to be passed to be added to each realisation, with one value per fault. "
        "The first argument should be the name of the value, "
        "the second the filepath to the space separated file containing the values. "
        "The file should have two columns, the name of a station followed by the value for that station, "
        "seperated by some number of spaces. "
        "If multiple source parameters are required this argument should be repeated.",
        default=[],
    )

    args = parser.parse_args()

    errors = []

    if not isfile(args.fault_selection_file):
        errors.append(f"Specified station file not found: {args.fault_selection_file}")
    if not isfile(args.gcmt_file):
        errors.append(f"Specified gcmt file not found: {args.gcmt_file}")

    if args.aggregate_file is not None and isfile(args.aggregate_file):
        errors.append(
            f"Specified aggregation file {args.aggregate_file} already exists, please choose another file"
        )
    if args.vel_mod_1d is not None and not isfile(args.vel_mod_1d):
        errors.append(
            f"Specified 1d velocity model file {args.vel_mod_1d} does not exist"
        )

    for i, (param_name, filepath) in enumerate(args.source_parameter):
        filepath = abspath(filepath)
        if not isfile(filepath):
            errors.append(
                f"The file {filepath} given for parameter "
                f"{param_name} does not exist"
            )
        else:
            args.source_parameter[i][1] = filepath

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


def generate_fault_realisations(
    data: GCMT_Source,
    realisation_count: int,
    cybershake_root: str,
    perturbation_function: Callable,
    aggregate_file: Union[str, None],
    vel_mod_1d: pd.DataFrame,
    primary_logger_name: str,
    additional_source_parameters: Dict[str, Any],
):
    primary_logger = get_logger(name=primary_logger_name)
    fault_logger = get_realisation_logger(primary_logger, data.pid)
    fault_logger.debug(f"Fault {data.pid} had data {data}")
    fault_name = data.pid

    for i in range(1, realisation_count + 1):
        realisation_name = get_realisation_name(fault_name, i)

        realisation_file_name = get_srf_path(cybershake_root, realisation_name)
        vel_mod_1d_dir = get_realisation_VM_dir(cybershake_root, realisation_name)
        makedirs(vel_mod_1d_dir, exist_ok=True)

        fault_logger.debug(
            f"Generating realisation {i} of {realisation_count} for fault {fault_name}"
        )

        generate_realisation(
            realisation_file_name,
            realisation_name,
            perturbation_function,
            data,
            additional_source_parameters,
            aggregate_file,
            vel_mod_1d,
            vel_mod_1d_dir,
            fault_logger,
        )

    unperturbated_function = load_perturbation_function(UNPERTURBATED)
    if perturbation_function != unperturbated_function:
        unperturbated_realisation = unperturbated_function(
            sources_line=data, additional_source_parameters=additional_source_parameters
        )
        rel_df = pd.DataFrame(unperturbated_realisation, index=[0])
        realisation_file_name = join(
            get_sources_dir(cybershake_root), fault_name, f"{fault_name}.csv"
        )
        rel_df.to_csv(realisation_file_name, index=False)


def generate_messages(
    additional_source_parameters: pd.DataFrame,
    aggregate_file,
    cybershake_root,
    faults,
    gcmt_lines,
    perturbation_function,
    vel_mod_1d,
    primary_logger,
):
    messages = []
    for i in range(len(faults)):
        gcmt_data = gcmt_lines[i]

        fault_name = gcmt_data[0]
        additional_source_specific_data = {}

        if fault_name in additional_source_parameters.index:
            additional_source_specific_data = additional_source_parameters.loc[
                fault_name
            ].to_dict()

        if vel_mod_1d is not None:
            vel_mod_1d_layers = load_1d_velocity_mod(vel_mod_1d)
        else:
            vel_mod_1d_layers = None

        messages.append(
            (
                GCMT_Source(*gcmt_data),
                faults[fault_name],
                cybershake_root,
                perturbation_function,
                aggregate_file,
                vel_mod_1d_layers,
                primary_logger.name,
                additional_source_specific_data,
            )
        )
    return messages


def main():

    primary_logger = get_logger("realisations_from_gcmt")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    faults = load_fault_selection_file(args.fault_selection_file)

    gcmt_data = pd.read_csv(
        args.gcmt_file, usecols=GCMT_FILE_COLUMNS, dtype=GCMT_DTYPE
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

    additional_source_parameters = pd.DataFrame()

    for param_name, filepath in args.source_parameter:
        parameter_df = pd.read_csv(
            filepath,
            delim_whitespace=True,
            header=None,
            index_col=0,
            names=[param_name],
            dtype={0: str},
        )
        additional_source_parameters = additional_source_parameters.join(
            parameter_df, how="outer"
        )

    messages = generate_messages(
        additional_source_parameters,
        args.aggregate_file,
        args.cybershake_root,
        faults,
        gcmt_lines,
        perturbation_function,
        args.vel_mod_1d,
        primary_logger,
    )

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
    worker_pool.starmap(generate_fault_realisations, messages)


if __name__ == "__main__":
    main()
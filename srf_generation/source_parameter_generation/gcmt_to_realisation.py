#! /usr/bin/python3

import argparse
from logging import Logger
from os import makedirs
from os.path import abspath, isfile, dirname, join
from typing import Callable, Union, Dict, Any, Tuple, List

import pandas as pd
from numpy import digitize

from qcore.simulation_structure import get_realisation_name
from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_realisation_logger,
)

from srf_generation.source_parameter_generation.uncertainties.common import (
    GCMT_PARAM_NAMES,
    GCMT_Source,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
)

DEFAULT_1D_VELOCITY_MODEL_PATH = join(
    dirname(__file__), "velocity_model", "lp_generic1d-gp01_v1.vmod"
)

GCMT_FILE_COLUMNS = [
    "PublicID",
    "Latitude",
    "Longitude",
    "strike1",
    "dip1",
    "rake1",
    "Mw",
    "CD",
]

UNPERTURBATED = "unperturbated"
GCMT_DTYPE = {
    GCMT_FILE_COLUMNS[i]: col_type
    for i, col_type in enumerate(
        [
            "object",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
            "float64",
        ]
    )
}


def load_args(primary_logger: Logger):
    parser = argparse.ArgumentParser(
        "Generates realisation files from gcmt and other input data."
    )

    parser.add_argument("fault_name")
    parser.add_argument("realisation_count", type=int)
    parser.add_argument("gcmt_file", type=abspath)
    parser.add_argument("--version", type=str, default="unperturbated")
    parser.add_argument(
        "--vel_mod_1d", type=abspath, default=DEFAULT_1D_VELOCITY_MODEL_PATH
    )
    parser.add_argument("--vel_mod_1d_out", type=abspath)
    parser.add_argument("--output_dir", "-o", type=abspath, default=abspath("."))
    parser.add_argument(
        "--aggregate_file",
        "-a",
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

    if not isfile(args.gcmt_file):
        errors.append(f"Specified gcmt file not found: {args.gcmt_file}")

    if args.aggregate_file is not None and isfile(args.aggregate_file):
        errors.append(
            f"Specified aggregation file {args.aggregate_file} already exists, please choose another file"
        )
    if not isfile(args.vel_mod_1d):
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

    makedirs(args.output_dir, exist_ok=True)

    log_file = join(args.output_dir, primary_logger.name + "_log.txt")
    add_general_file_handler(primary_logger, log_file)

    primary_logger.debug(f"Arguments parsed: {args}")

    return args


def load_1d_velocity_mod(vel_mod_1d):
    return pd.read_csv(
        vel_mod_1d,
        delim_whitespace=True,
        header=None,
        skiprows=1,
        names=["depth", "vp", "vs", "rho", "qp", "qs"],
    )


def save_1d_velocity_model(
    perturbed_realisation: pd.DataFrame, vel_mod_1d_dir, realisation_name
):
    file_name = join(vel_mod_1d_dir, f"{realisation_name}.v1d")
    with open(file_name, "w") as v1d:
        v1d.write(f"{len(perturbed_realisation)}\n")
        for line in perturbed_realisation.itertuples():
            v1d.write(
                f"{line.depth} {line.vp} {line.vs} {line.rho} {line.qp} {line.qs}\n"
            )
    return file_name


def generate_fault_realisations(
    data: GCMT_Source,
    realisation_count: int,
    output_directory: str,
    perturbation_function: Callable,
    aggregate_file: Union[str, None],
    vel_mod_1d: pd.DataFrame,
    vel_mod_1d_dir: str,
    primary_logger_name: str,
    additional_source_parameters: Dict[str, Any],
):
    primary_logger = get_logger(name=primary_logger_name)
    fault_logger = get_realisation_logger(primary_logger, data.pid)
    fault_logger.debug(f"Fault {data.pid} had data {data}")
    fault_name = data.pid

    for i in range(1, realisation_count + 1):
        realisation_name = get_realisation_name(fault_name, i)
        realisation_file_name = join(output_directory, f"{realisation_name}.csv")
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
            sources_line=data,
            additional_source_parameters=additional_source_parameters,
            vel_mod_1d=None,
        )
        rel_df = pd.DataFrame(unperturbated_realisation, index=[0])
        realisation_file_name = join(output_directory, f"{fault_name}.csv")
        rel_df.to_csv(realisation_file_name, index=False)


def generate_realisation(
    realisation_file_name,
    realisation_name,
    perturbation_function,
    data,
    additional_source_parameters,
    aggregate_file,
    vel_mod_1d,
    vel_mod_1d_dir,
    fault_logger,
):
    fault_logger.debug("Calling perturbation function.")
    perturbed_realisation = perturbation_function(
        sources_line=data,
        additional_source_parameters=additional_source_parameters,
        vel_mod_1d=vel_mod_1d,
    )
    fault_logger.debug(
        f"Got results from perturbation_function: {perturbed_realisation}"
    )
    perturbed_realisation["params"]["name"] = realisation_name

    if (
        vel_mod_1d_dir is not None
        and "vel_mod_1d" in perturbed_realisation.keys()
        and not vel_mod_1d.equals(perturbed_realisation["vel_mod_1d"])
    ):
        perturbed_vel_mod_1d = perturbed_realisation["vel_mod_1d"]
        makedirs(vel_mod_1d_dir, exist_ok=True)
        file_name_1d_vel_mod = save_1d_velocity_model(
            perturbed_vel_mod_1d, vel_mod_1d_dir, realisation_name
        )
        perturbed_realisation["params"]["v_mod_1d_name"] = file_name_1d_vel_mod

    makedirs(dirname(realisation_file_name), exist_ok=True)
    fault_logger.debug(
        f"Created Srf directory and attempting to save perturbated source generation parameters there: {realisation_file_name}"
    )
    rel_df = pd.DataFrame(perturbed_realisation["params"], index=[0])
    rel_df.to_csv(realisation_file_name, index=False)

    if aggregate_file is not None:
        if not isfile(aggregate_file):
            rel_df.to_csv(aggregate_file)
        else:
            rel_df.to_csv(aggregate_file, mode="a", header=False)

    fault_logger.debug(
        f"Parameters saved succesfully. Continuing to next realisation if one exists."
    )


def get_additional_source_parameters(
    source_parameters: List[Tuple[str, str]],
    gcmt_data: pd.DataFrame,
    vel_mod_1d_layers: pd.DataFrame,
):
    """Gets the sheer wave velocities at the point source depth and
    also takes in the source parameter name-filepath pairs passed as arguments,
    extracting the values for each event
    :param source_parameters: A list of tuples containing name, filepath pairs
    :param gcmt_data:
    :param vel_mod_1d_layers: A dataframe containing a row for every layer of the velocity model
    """
    # For each event get the layer it corresponds to by taking the cumulative sum of the layer depths and finding
    # the layer with bottom depth lower than the event
    depth_bins = digitize(
        gcmt_data["depth"].round(5), vel_mod_1d_layers["depth"].cumsum().round(5)
    )
    additional_source_parameters = pd.DataFrame(
        {"vs": vel_mod_1d_layers["vs"].iloc[depth_bins].values}, gcmt_data["pid"].values
    )
    for param_name, filepath in source_parameters:
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
    return additional_source_parameters


def main():

    primary_logger = get_logger("GCMT_2_realisation")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    gcmt_data = pd.read_csv(args.gcmt_file, usecols=GCMT_FILE_COLUMNS, index_col=0)[
        ["Latitude", "Longitude", "CD", "Mw", "strike1", "dip1", "rake1"]
    ]
    gcmt_data.columns = GCMT_PARAM_NAMES[1:]

    primary_logger.debug(
        f"GCMT data for {len(gcmt_data)} faults were loaded. Checking for fault {args.fault_name}"
    )

    if args.fault_name in gcmt_data.index:
        gcmt_line = gcmt_data.loc[args.fault_name]
    else:
        message = (
            f"Fault {args.fault_name} is not available in the given GCMT file, exiting"
        )
        primary_logger.log(NOPRINTCRITICAL, message)
        primary_logger.debug("Available faults are as follows:")
        primary_logger.debug(gcmt_data.index)
        raise ValueError(message)

    vel_mod_1d_layers = load_1d_velocity_mod(args.vel_mod_1d)

    additional_source_parameters = get_additional_source_parameters(
        args.source_parameter, gcmt_data, vel_mod_1d_layers
    )

    if gcmt_line[0] in additional_source_parameters.index:
        additional_source_specific_data = additional_source_parameters.loc[
            gcmt_line[0]
        ].to_dict()
    else:
        additional_source_specific_data = {}

    vel_mod_1d_layers = load_1d_velocity_mod(args.vel_mod_1d)

    if args.realisation_count > 1:
        generate_fault_realisations(
            GCMT_Source(args.fault_name, *gcmt_line),
            args.realisation_count,
            args.output_dir,
            perturbation_function,
            args.aggregate_file,
            vel_mod_1d_layers,
            args.vel_mod_1d_out,
            primary_logger.name,
            additional_source_specific_data,
        )
    else:
        generate_realisation(
            join(args.output_dir, f"{args.fault_name}.csv"),
            args.fault_name,
            perturbation_function,
            GCMT_Source(args.fault_name, *gcmt_line),
            additional_source_parameters,
            args.aggregate_file,
            args.vel_mod_1d,
            args.vel_mod_1d_out,
            primary_logger,
        )


if __name__ == "__main__":
    main()

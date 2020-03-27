#! /usr/bin/python3

import argparse
from logging import Logger
from os import makedirs
from os.path import abspath, isfile, dirname, join
from typing import Callable, Union, Dict, Any, Tuple, List

import pandas as pd
from qcore.nhm import load_nhm, NHMFault
from qcore.formats import load_vs30_file

from qcore.simulation_structure import get_realisation_name
from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_realisation_logger,
    get_basic_logger,
)

from srf_generation.source_parameter_generation.uncertainties.common import (
    NHM_Source,
    get_seed,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
)

DEFAULT_1D_VELOCITY_MODEL_PATH = join(
    dirname(__file__), "velocity_model", "lp_generic1d-gp01_v1.vmod"
)


def load_args(primary_logger: Logger):
    parser = argparse.ArgumentParser(
        "Generates realisation files from nhm and other input data."
    )

    parser.add_argument("fault_name")
    parser.add_argument("realisation_count", type=int)
    parser.add_argument(
        "nhm_file",
        type=abspath,
        help="The path to a national seismic hazard model file. "
        "Must have entries for all events named in the fault selection file. "
        "Additional events not named will be ignored.",
    )
    parser.add_argument("type", type=str, help="The type of srf to generate.")

    add_common_arguments(parser)

    args = parser.parse_args()

    errors = []

    verify_args(args, errors, primary_logger)

    if errors:
        message = (
            "At least one error was detected when verifying arguments:\n"
            + "\n".join(errors)
        )
        primary_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)

    makedirs(args.output_dir, exist_ok=True)

    log_file = join(args.output_dir, primary_logger.name + "_log.txt")
    add_general_file_handler(primary_logger, log_file)

    primary_logger.debug(f"Arguments parsed: {args}")

    return args


def verify_args(args, errors, parser_logger=get_basic_logger()):

    if args.version is None:
        if args.type is not None:
            args.version = f"nhm_{args.type}"
        else:
            parser_logger.debug(
                "No version or type given, generating type 1 realisations"
            )
            args.version = f"nhm_4"

    if not isfile(args.nhm_file):
        errors.append(f"Specified nhm file not found: {args.nhm_file}")
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
    if args.vs30_median is not None:
        if not isfile(args.vs30_median):
            errors.append(
                f"The file {args.vs30_median} given for parameter --vs30_median does not exist"
            )
        if args.vs30_sigma is not None and not isfile(args.vs30_sigma):
            errors.append(
                f"The file {args.vs30_sigma} given for parameter --vs30_sigma does not exist"
            )
    elif args.vs30_sigma is not None:
        errors.append(
            f"If the vs30 sigma file is given the median file should also be given"
        )


def add_common_arguments(parser, single_event=True):
    parser.add_argument(
        "--version",
        type=str,
        help="The name of the version to perturbate the input parameters with. "
        "By default no perturbation will occur. "
        "Should be the name of the file without the .py suffix.",
    )
    parser.add_argument(
        "--vel_mod_1d", type=abspath, default=DEFAULT_1D_VELOCITY_MODEL_PATH
    )
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
        "separated by some number of spaces. "
        "If multiple source parameters are required this argument should be repeated.",
        default=[],
    )
    parser.add_argument(
        "--common_source_parameter",
        type=str,
        nargs=2,
        action="append",
        metavar=("name", "parameter"),
        help="Values to be passed to be added to each realisation, with the same value for every event. "
        "The first argument should be the name of the value, the second should be the value. "
        "If the value is a valid number it will be treated as a float, otherwise it will be a string"
        "If multiple source parameters are required this argument should be repeated.",
        default=[],
    )
    vs30_parser = parser.add_argument_group(
        title="vs30 arguments",
        description="All arguments in this group must be given if any are given",
    )
    vs30_parser.add_argument(
        "--vs30_median",
        type=abspath,
        help="The path to a file containing median VS30s.",
    )
    vs30_parser.add_argument(
        "--vs30_sigma", type=abspath, help="The path to a file containing VS30 sigmas."
    )
    if single_event:
        parser.add_argument("--output_dir", "-o", type=abspath, default=abspath("."))
        parser.add_argument("--vel_mod_1d_out", type=abspath)
        vs30_parser.add_argument(
            "--vs30_out",
            type=abspath,
            help="The path to a file to save the perturbated VS30s to",
            default=abspath("."),
        )
    else:
        parser.add_argument(
            "--cybershake_root",
            type=abspath,
            default=abspath("."),
            help="The path to the root of the simulation root directory. Defaults to the current directory.",
        )


def load_vs30_median_sigma(vs30_median, vs30_sigma):
    if vs30_median is None:
        return None
    vs30_medians = load_vs30_file(vs30_median)
    vs30_medians.columns = ["median"]
    if vs30_sigma is not None:
        vs30_sigmas: pd.DataFrame = load_vs30_file(vs30_sigma)
        vs30_sigmas.columns = ["sigma"]
        vs30 = vs30_medians.join(vs30_sigmas)
    else:
        vs30 = vs30_medians
        vs30["sigma"] = 0.0
    return vs30


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
    data: NHMFault,
    realisation_count: int,
    output_directory: str,
    perturbation_function: Callable,
    unperturbation_function: Callable,
    aggregate_file: Union[str, None],
    vel_mod_1d: pd.DataFrame,
    vel_mod_1d_dir: str,
    vs30_data: pd.DataFrame,
    vs30_out_file: str,
    primary_logger_name: str,
    additional_source_parameters: Dict[str, Any],
):
    primary_logger = get_logger(name=primary_logger_name)
    fault_logger = get_realisation_logger(primary_logger, data.name)
    fault_logger.debug(f"Fault {data.name} had data {data}")
    fault_name = data.name

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
            vs30_data,
            vs30_out_file,
            fault_logger,
        )

    if perturbation_function != unperturbation_function:
        unperturbated_realisation = unperturbation_function(
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
    data: NHMFault,
    additional_source_parameters,
    aggregate_file,
    vel_mod_1d,
    vel_mod_1d_dir,
    vs30_data: pd.DataFrame,
    vs30_out_file: str,
    fault_logger: Logger = get_basic_logger(),
):
    fault_logger.debug("Calling perturbation function.")
    perturbed_realisation = perturbation_function(
        source_data=data,
        additional_source_parameters=additional_source_parameters,
        vel_mod_1d=vel_mod_1d,
        vs30_data=vs30_data,
    )
    fault_logger.debug(
        f"Got results from perturbation_function: {perturbed_realisation}"
    )
    perturbed_realisation["params"]["name"] = realisation_name
    if "srfgen_seed" not in perturbed_realisation["params"]:
        perturbed_realisation["params"]["srfgen_seed"] = get_seed()

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

    if vs30_out_file is not None and "vs30" in perturbed_realisation.keys():
        perturbated_vs30: pd.DataFrame = perturbed_realisation["vs30"]
        perturbated_vs30.to_csv(
            vs30_out_file, columns="vs30", sep=" ", index=True, header=False
        )
        perturbed_realisation["params"]["vs30_file_path"] = vs30_out_file

    if "z_values" in perturbed_realisation.keys():
        z_df = pd.DataFrame(perturbed_realisation["z_values"], index=[0])
        z_df.to_csv(realisation_file_name.replace(".csv", "_z_values.csv"), index=False)

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
        f"Parameters saved succesfully to {realisation_file_name}. Continuing to next realisation if one exists."
    )


def get_additional_source_parameters(
    source_parameters: List[Tuple[str, str]],
    common_source_parameters: List[Tuple[str, str]],
    nhm_data: Dict[str, NHMFault],
    vel_mod_1d_layers: pd.DataFrame,
):
    """Takes in the source parameter name-filepath pairs passed as arguments,
    extracting the values for each event
    :param source_parameters: A list of tuples containing name, filepath pairs
    :param common_source_parameters: A list of tuples containing name, value pairs
    :param nhm_data: A Dictionary of faults and their NHMFault objects
    :param vel_mod_1d_layers: A dataframe containing a row for every layer of the velocity model
    """
    additional_source_parameters = pd.DataFrame(index=sorted(list(nhm_data.keys())),)
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
    for param_name, value in common_source_parameters:
        try:
            value = int(value)
        except:
            try:
                value = float(value)
            except:
                pass
        additional_source_parameters[param_name] = value
    return additional_source_parameters


def main():

    primary_logger = get_logger("NHM_2_realisation")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    unperturbation_function = load_perturbation_function(f"nhm_{args.type}")
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    nhm_data = load_nhm(args.nhm_file)

    primary_logger.debug(
        f"NHM data for {len(nhm_data.keys())} faults were loaded. Checking for fault {args.fault_name}"
    )

    if args.fault_name in nhm_data.keys():
        fault_nhm = nhm_data[args.fault_name]
    else:
        message = (
            f"Fault {args.fault_name} is not available in the given NHM file, exiting"
        )
        primary_logger.log(NOPRINTCRITICAL, message)
        primary_logger.debug("Available faults are as follows:")
        primary_logger.debug(nhm_data.keys())
        raise ValueError(message)

    vel_mod_1d_layers = load_1d_velocity_mod(args.vel_mod_1d)

    additional_source_parameters = get_additional_source_parameters(
        args.source_parameter,
        args.common_source_parameter,
        fault_nhm,
        vel_mod_1d_layers,
    )

    vs30 = load_vs30_median_sigma(args.vs30_median, args.vs30_sigma)

    if fault_nhm.name in additional_source_parameters.index:
        additional_source_specific_data = additional_source_parameters.loc[
            fault_nhm.name
        ].to_dict()
    else:
        additional_source_specific_data = {}

    vel_mod_1d_layers = load_1d_velocity_mod(args.vel_mod_1d)

    if args.realisation_count > 1:
        generate_fault_realisations(
            fault_nhm,
            args.realisation_count,
            args.output_dir,
            perturbation_function,
            unperturbation_function,
            args.aggregate_file,
            vel_mod_1d_layers,
            args.vel_mod_1d_out,
            vs30,
            args.vs30_out,
            primary_logger.name,
            additional_source_specific_data,
        )
    else:
        generate_realisation(
            join(args.output_dir, f"{args.fault_name}.csv"),
            args.fault_name,
            perturbation_function,
            fault_nhm,
            additional_source_parameters,
            args.aggregate_file,
            args.vel_mod_1d,
            args.vel_mod_1d_out,
            vs30,
            args.vs30_out,
            primary_logger,
        )


if __name__ == "__main__":
    main()

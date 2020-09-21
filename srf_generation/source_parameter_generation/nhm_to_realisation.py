#! /usr/bin/python3

import argparse
from logging import Logger
from os import makedirs
from os.path import abspath, isfile, join
from typing import Callable, Union, Dict, Any, Tuple, List

import pandas as pd
from qcore.nhm import load_nhm, NHMFault

from qcore.simulation_structure import get_realisation_name
from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_realisation_logger,
    get_basic_logger,
)

from srf_generation.source_parameter_generation.common import (
    add_common_arguments,
    load_vs30_median_sigma,
    load_1d_velocity_mod,
    process_perturbation,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
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
                "No version or type given, generating type 4 realisations"
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
            None,
            vs30_data,
            vs30_out_file,
            fault_logger,
        )

    if perturbation_function != unperturbation_function:
        unperturbated_realisation = unperturbation_function(
            source_data=data,
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
    hf_vel_mod_1d,
    vs30_data: pd.DataFrame,
    vs30_out_file: str,
    fault_logger: Logger = get_basic_logger(),
):
    fault_logger.debug("Calling perturbation function.")
    perturbed_realisation = perturbation_function(
        source_data=data,
        additional_source_parameters=additional_source_parameters.copy(),
        vel_mod_1d=vel_mod_1d.copy(deep=True) if vel_mod_1d is not None else None,
        vs30_data=vs30_data.copy(deep=True) if vs30_data is not None else None,
    )
    fault_logger.debug(
        f"Got results from perturbation_function: {perturbed_realisation}"
    )
    process_perturbation(
        aggregate_file,
        fault_logger,
        hf_vel_mod_1d,
        perturbed_realisation,
        realisation_file_name,
        realisation_name,
        vel_mod_1d,
        vel_mod_1d_dir,
        vs30_out_file,
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
    additional_source_parameters = pd.DataFrame(
        index=sorted(list(nhm_data.keys())),
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
            None,
            vs30,
            args.vs30_out,
            primary_logger,
        )


if __name__ == "__main__":
    main()

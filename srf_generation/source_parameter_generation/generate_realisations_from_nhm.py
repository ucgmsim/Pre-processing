#! /usr/bin/python3

import argparse
from logging import Logger
from multiprocessing import pool
from os import makedirs
from os.path import abspath, isfile, join
from typing import Callable, Union, Dict, Any

import pandas as pd

from qcore.formats import load_fault_selection_file
from qcore.nhm import NHMFault, load_nhm
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
from srf_generation.source_parameter_generation.nhm_to_realisation import (
    generate_realisation,
    get_additional_source_parameters,
    verify_args,
)
from srf_generation.source_parameter_generation.common import (
    add_common_arguments,
    load_vs30_median_sigma,
    load_1d_velocity_mod,
)

from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
)


def load_args(primary_logger: Logger):
    parser = argparse.ArgumentParser(
        "Generates realisation files from nhm and other input data. "
        "Intended to be used with the standard qcore simultaion structure."
    )

    parser.add_argument(
        "fault_selection_file",
        type=abspath,
        help="The path to a white space delimited two column file containing a list of events, "
        "followed by the number of realisations for that event. "
        "If any events are missing from the nhm file the program will stop. "
        "If any events have only one realisation they will not have the _RELXX suffix attached.",
    )
    parser.add_argument(
        "nhm_file",
        type=abspath,
        help="The path to a national seismic hazard model file. "
        "Must have entries for all events named in the fault selection file. "
        "Additional events not named will be ignored.",
    )
    parser.add_argument("type", type=str, help="The type of srf to generate.")
    parser.add_argument(
        "--n_processes",
        type=int,
        help="The number of processes to run at once. Capped at the number of events to generate realisations for.",
        default=1,
    )

    add_common_arguments(parser, single_event=False)

    args = parser.parse_args()
    primary_logger.debug(f"Raw arguments passed, beginning argument processing: {args}")

    errors = []

    if not isfile(args.fault_selection_file):
        errors.append(
            f"Specified selection file not found: {args.fault_selection_file}"
        )

    verify_args(args, errors, primary_logger)

    if errors:
        message = (
            "At least one error was detected when verifying arguments:\n"
            + "\n".join(errors)
        )
        primary_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)

    makedirs(args.cybershake_root, exist_ok=True)

    log_file = join(args.cybershake_root, primary_logger.name + "_log.txt")
    add_general_file_handler(primary_logger, log_file)

    primary_logger.debug(f"Arguments parsed: {args}")

    return args


def generate_fault_realisations(
    data: NHMFault,
    realisation_count: int,
    cybershake_root: str,
    perturbation_function: Callable,
    unperturbation_function: Callable,
    aggregate_file: Union[str, None],
    vel_mod_1d: pd.DataFrame,
    vs30_data: pd.DataFrame,
    primary_logger_name: str,
    additional_source_parameters: Dict[str, Any],
):
    primary_logger = get_logger(name=primary_logger_name)
    fault_logger = get_realisation_logger(primary_logger, data.name)
    fault_logger.debug(f"Fault {data.name} had data {data}")
    fault_name = data.name
    fault_logger.info(f"Generating realisations for event {fault_name}")

    if realisation_count == 1:
        fault_logger.debug(f"Generating the only realisation of fault {fault_name}")
        vs30_out_file = join(
            get_realisation_VM_dir(cybershake_root, fault_name), f"{fault_name}.vs30"
        )
        generate_realisation(
            get_srf_path(cybershake_root, fault_name).replace(".srf", ".csv"),
            fault_name,
            perturbation_function,
            data,
            additional_source_parameters,
            aggregate_file,
            vel_mod_1d,
            get_realisation_VM_dir(cybershake_root, fault_name),
            vs30_data,
            vs30_out_file,
            fault_logger,
        )
        return

    for i in range(1, realisation_count + 1):
        realisation_name = get_realisation_name(fault_name, i)

        realisation_file_name = get_srf_path(cybershake_root, realisation_name).replace(
            ".srf", ".csv"
        )

        # if isfile(realisation_file_name):
        #     continue

        vel_mod_1d_dir = get_realisation_VM_dir(cybershake_root, realisation_name)
        vs30_out_file = join(vel_mod_1d_dir, f"{realisation_name}.vs30")

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
            source_data=data,
            additional_source_parameters=additional_source_parameters,
            vel_mod_1d=None,
        )
        rel_df = pd.DataFrame(unperturbated_realisation["params"], index=[0])
        realisation_file_name = join(
            get_sources_dir(cybershake_root), fault_name, f"{fault_name}.csv"
        )
        rel_df.to_csv(realisation_file_name, index=False)


def generate_messages(
    additional_source_parameters: pd.DataFrame,
    aggregate_file,
    cybershake_root,
    faults,
    nhm_faults,
    perturbation_function,
    unperturbation_function,
    vel_mod_1d,
    vs30_data: pd.DataFrame,
    primary_logger,
):
    messages = []
    for fault_name, realisation_count in faults.items():
        fault_data = nhm_faults[fault_name]

        additional_source_specific_data = {}

        if fault_name in additional_source_parameters.index:
            additional_source_specific_data = additional_source_parameters.loc[
                fault_name
            ].to_dict()

        messages.append(
            (
                fault_data,
                faults[fault_name],
                cybershake_root,
                perturbation_function,
                unperturbation_function,
                aggregate_file,
                vel_mod_1d,
                vs30_data,
                primary_logger.name,
                additional_source_specific_data,
            )
        )
    return messages


def main():

    primary_logger = get_logger("realisations_from_nhm")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    unperturbation_function = load_perturbation_function(f"nhm_{args.type}")
    primary_logger.debug(f"Perturbation function loaded. Version: {args.version}")

    faults = load_fault_selection_file(args.fault_selection_file)

    nhm_data = load_nhm(args.nhm_file)

    primary_logger.debug(
        f"{len(faults)} faults were selected and nhm data for {len(nhm_data)} faults were loaded. "
        f"Running cross comparison now."
    )

    nhm_faults = {key: nhm_data.get(key) for key in faults.keys()}

    missing_faults = [key for key in nhm_faults.keys() if nhm_faults[key] is None]

    if missing_faults:
        message = f"Some fault(s) in the fault selection were not found in the nhm file: {missing_faults}"
        primary_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)

    primary_logger.debug(
        f"All {len(nhm_faults.keys())} faults specified in the fault selection file were found"
    )

    velocity_model_1d = load_1d_velocity_mod(args.vel_mod_1d)
    vs30 = load_vs30_median_sigma(args.vs30_median, args.vs30_sigma)

    additional_source_parameters = get_additional_source_parameters(
        args.source_parameter,
        args.common_source_parameter,
        nhm_faults,
        velocity_model_1d,
    )

    messages = generate_messages(
        additional_source_parameters,
        args.aggregate_file,
        args.cybershake_root,
        faults,
        nhm_faults,
        perturbation_function,
        unperturbation_function,
        velocity_model_1d,
        vs30,
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

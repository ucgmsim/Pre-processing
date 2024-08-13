#! /usr/bin/python3

import argparse
from logging import Logger
from os import makedirs
from os.path import abspath, isfile, dirname, join
from typing import Tuple, List

import pandas as pd

from qcore.qclogging import (
    get_logger,
    add_general_file_handler,
    NOPRINTCRITICAL,
    get_basic_logger,
)

from srf_generation.source_parameter_generation.common import (
    save_1d_velocity_model,
    load_vs30_median_sigma,
    load_1d_velocity_mod,
    add_common_arguments,
    get_depth_property,
    write_asperites,
    generate_fault_realisations,
)
from srf_generation.source_parameter_generation.uncertainties.common import (
    GCMT_PARAM_NAMES,
    GCMT_Source,
    get_seed,
)
from srf_generation.source_parameter_generation.uncertainties.versions import (
    load_perturbation_function,
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
    parser.add_argument(
        "gcmt_file",
        type=abspath,
        help="The path to a geonet cmt solutions file. "
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
            args.version = f"gcmt_{args.type}"
        else:
            parser_logger.debug(
                "No version or type given, generating type 1 realisations"
            )
            args.version = f"gcmt_1"

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


def generate_realisation(
    realisation_file_name,
    realisation_name,
    perturbation_function,
    data,
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
        perturbed_vel_mod_1d = perturbed_realisation.pop("vel_mod_1d")
        makedirs(vel_mod_1d_dir, exist_ok=True)
        file_name_srf_1d_vel_mod = save_1d_velocity_model(
            perturbed_vel_mod_1d, vel_mod_1d_dir, f"srf_{realisation_name}"
        )
        perturbed_realisation["params"]["srf_vel_mod_1d"] = file_name_srf_1d_vel_mod

    if (
        vel_mod_1d_dir is not None
        and "hf_vel_mod_1d" in perturbed_realisation.keys()
        and (
            hf_vel_mod_1d is None
            or not hf_vel_mod_1d.equals(perturbed_realisation["hf_vel_mod_1d"])
        )
    ):
        perturbed_vel_mod_1d = perturbed_realisation.pop("hf_vel_mod_1d")
        makedirs(vel_mod_1d_dir, exist_ok=True)
        file_name_1d_vel_mod = save_1d_velocity_model(
            perturbed_vel_mod_1d, vel_mod_1d_dir, realisation_name
        )
        perturbed_realisation["params"]["v_mod_1d_name"] = file_name_1d_vel_mod

    if vs30_out_file is not None and "vs30" in perturbed_realisation.keys():
        perturbed_vs30: pd.DataFrame = perturbed_realisation.pop("vs30")
        makedirs(dirname(vs30_out_file), exist_ok=True)
        perturbed_vs30.to_csv(
            vs30_out_file, columns=["vs30"], sep=" ", index=True, header=False
        )
        perturbed_realisation["params"]["vs30_file_path"] = vs30_out_file

    if "z_values" in perturbed_realisation.keys():
        z_df = pd.DataFrame(perturbed_realisation["z_values"], index=[0])
        z_df.to_csv(realisation_file_name.replace(".csv", "_z_values.csv"), index=False)

    makedirs(dirname(realisation_file_name), exist_ok=True)

    if "asperities" in perturbed_realisation.keys():
        asperity_file_path = realisation_file_name.replace(".csv", ".aspf")
        write_asperites(perturbed_realisation["asperities"], asperity_file_path)
        perturbed_realisation["params"]["asperity_file"] = asperity_file_path

    fault_logger.debug(
        f"Created Srf directory and attempting to save perturbed source generation parameters there: {realisation_file_name}"
    )
    rel_df = pd.DataFrame(perturbed_realisation["params"], index=[0])
    rel_df.to_csv(realisation_file_name, index=False)

    if aggregate_file is not None:
        if not isfile(aggregate_file):
            rel_df.to_csv(aggregate_file, index=False)
        else:
            rel_df.to_csv(aggregate_file, mode="a", header=False, index=False)

    fault_logger.debug(
        f"Parameters saved successfully. Continuing to next realisation if one exists."
    )


def get_additional_source_parameters(
    source_parameters: List[Tuple[str, str]],
    common_source_parameters: List[Tuple[str, str]],
    gcmt_data: pd.DataFrame,
    vel_mod_1d_layers: pd.DataFrame,
):
    """Gets the shear wave velocities at the point source depth and
    also takes in the source parameter name-filepath pairs passed as arguments,
    extracting the values for each event
    :param source_parameters: A list of tuples containing name, filepath pairs
    :param common_source_parameters: A list of tuples containing name, value pairs
    :param gcmt_data:
    :param vel_mod_1d_layers: A dataframe containing a row for every layer of the velocity model
    """
    # For each event get the layer it corresponds to by taking the cumulative sum of the layer depths and finding
    # the layer with bottom depth lower than the event
    additional_source_parameters = pd.DataFrame(
        {
            "vs": get_depth_property(gcmt_data["depth"], vel_mod_1d_layers, "vs"),
            "rho": get_depth_property(gcmt_data["depth"], vel_mod_1d_layers, "rho"),
        },
        gcmt_data["pid"].values,
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
    primary_logger = get_logger("GCMT_2_realisation")

    args = load_args(primary_logger)

    perturbation_function = load_perturbation_function(args.version)
    unperturbed_function = load_perturbation_function(f"gcmt_{args.version}")
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
        args.source_parameter,
        args.common_source_parameter,
        gcmt_data,
        vel_mod_1d_layers,
    )

    vs30 = load_vs30_median_sigma(args.vs30_median, args.vs30_sigma)

    if gcmt_line[0] in additional_source_parameters.index:
        additional_source_specific_data = additional_source_parameters.loc[
            gcmt_line[0]
        ].to_dict()
    else:
        additional_source_specific_data = {}

    vel_mod_1d_layers = load_1d_velocity_mod(args.vel_mod_1d)

    generate_fault_realisations(
        args.fault_name,
        GCMT_Source(args.fault_name, *gcmt_line),
        args.realisation_count,
        args.output_dir,
        perturbation_function,
        unperturbed_function,
        args.aggregate_file,
        vel_mod_1d_layers,
        args.vel_mod_1d_out,
        vs30,
        args.vs30_out,
        primary_logger.name,
        additional_source_specific_data,
        generate_realisation,
    )


if __name__ == "__main__":
    main()

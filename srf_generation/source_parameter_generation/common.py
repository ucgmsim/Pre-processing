from os.path import abspath, join, dirname
from typing import Callable, Union, Dict, Any

import pandas as pd
import numpy as np
from qcore.formats import load_vs30_file
from qcore.qclogging import get_logger, get_realisation_logger
from qcore.simulation_structure import get_realisation_name


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
        parser.add_argument(
            "-c",
            "--checkpointing",
            action="store_true",
            help="Activate checkpointing for realisation file generation. Any realisation files that already exist will not be regenerated. Default is to regenerate.",
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


DEFAULT_1D_VELOCITY_MODEL_PATH = join(
    dirname(__file__), "velocity_model", "lp_generic1d-gp01_v1.vmod"
)


def get_depth_property(target_depths, vel_mod_1d_layers, property_name):
    binned_depths = np.digitize(
        np.asarray(target_depths).round(5), vel_mod_1d_layers["depth"].cumsum().round(5)
    )
    return np.asarray(vel_mod_1d_layers[property_name].iloc[binned_depths])


def write_asperites(asperities_dict, output_file):
    background_value = asperities_dict["background"]
    asperities_list = asperities_dict["asperities"]
    with open(output_file, "w") as aspf:
        aspf.write(f"{background_value}\n")
        for asperity in asperities_list:
            aspf.write(f"{asperity.to_genslip_format()}\n")


def generate_fault_realisations(
    fault_name: str,
    data: Union["GCMT_Source", "NHM_data"],
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
    generate_realisation: Callable,
):
    primary_logger = get_logger(name=primary_logger_name)
    fault_logger = get_realisation_logger(primary_logger, fault_name)
    fault_logger.debug(f"Fault {fault_name} had data {data}")

    generate_realisation(
        join(output_directory, f"{fault_name}.csv"),
        fault_name,
        unperturbation_function,
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
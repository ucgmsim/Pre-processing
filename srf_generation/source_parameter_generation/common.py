from os import makedirs
from os.path import abspath, join, dirname, isfile

import pandas as pd
import numpy as np
from qcore.formats import load_vs30_file

from srf_generation.source_parameter_generation.uncertainties.common import get_seed


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


def process_perturbation(
    aggregate_file,
    fault_logger,
    hf_vel_mod_1d,
    perturbed_realisation,
    realisation_file_name,
    realisation_name,
    vel_mod_1d,
    vel_mod_1d_dir,
    vs30_out_file,
):
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
        perturbated_vs30: pd.DataFrame = perturbed_realisation.pop("vs30")
        perturbated_vs30.to_csv(
            vs30_out_file, columns=["vs30"], sep=" ", index=True, header=False
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

import argparse
from os.path import abspath
import pandas as pd

from qcore.vm_file import create_constant_vm_file
from qcore.utils import load_yaml

from srf_generation.velocity_model_generation.fault_damage_zone import (
    apply_fault_damage_zone,
    add_fault_damage_zone_properties,
)
from srf_generation.velocity_model_generation.generate_3d_velocity_model_perturbation import (
    generate_velocity_model_perturbation_file_from_config,
    load_parameter_file,
    generate_velocity_model_perturbation_file_from_model,
)


def load_args():
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Generates a single perturbation file from given inputs.",
    )
    parser.add_argument(
        "perturbation_file",
        type=abspath,
        help="The location to save the perturbation output file to.",
    )
    parser.add_argument(
        "vm_params_location",
        type=abspath,
        help="The location of the fault/events vm_params file.",
    )
    parser.add_argument(
        "-n",
        "--n_processes",
        default=1,
        type=int,
        help="The number of layers to generate simultaneously. Setting a value greater than the layer count will have no effect.",
    )
    parser.add_argument("-v", "--verbose", help="Output debug info to std_err")

    parser.add_argument(
        "--perturbation",
        action="store_true",
        help="Generate perturbated layers using Rob Graves Fractal3D code. If not used a base perturbation file of 1s will be used",
    )
    parser.add_argument(
        "--fault_damage_zone",
        action="store_true",
        help="Add a fault damage zone to the velocity model.",
    )

    perturbation_group = parser.add_argument_group(
        title="Perturbation",
        description="The settings to use for generating perturbations. Only one option may be used. File formats are available on the wiki.",
    )
    perturbation_group.add_argument(
        "--parameter_file",
        type=abspath,
        help="The location of a perturbation parameter file.",
    )
    perturbation_group.add_argument(
        "--model", type=abspath, help="The location of a perturbation model file."
    )

    fault_damage_zone_group = parser.add_argument_group(
        title="Fault damage zone", description="Parameters for a fault damage zone"
    )
    fault_damage_zone_group.add_argument(
        "--srf_location",
        type=abspath,
        help="Path to the srf. Required if --fault_damage_zone is passed, ignored otherwise",
    )
    add_fault_damage_zone_properties(fault_damage_zone_group)

    args = parser.parse_args()

    errors = []

    if args.perturbation:
        if (args.parameter_file is None) == (args.model is None):
            errors.append(
                "If --perturbation is given, one of --parameter_file and --model must also be given"
            )

    if args.fault_damage_zone:
        if args.srf_location is None:
            errors.append(
                "If --fault_damage_zone is given, --srf_location must also be given"
            )

    if len(errors) > 0:
        parser.error("\n".join(errors))

    return args


def main():
    """
    Combines the velocity model perturbation file functionality into one location
    """
    args = load_args()

    perturbation_file = args.perturbation_file
    vm_params = load_yaml(args.vm_params_location)
    processes = args.n_processes
    verbose = args.verbose

    if args.perturbation:
        if args.model:
            perturbation_model = pd.read_csv(args.model)
            generate_velocity_model_perturbation_file_from_model(
                vm_params, perturbation_model, perturbation_file, processes, verbose
            )
        elif args.parameter_file:
            common_params, layer_params = load_parameter_file(args.parameter_file)
            generate_velocity_model_perturbation_file_from_config(
                common_params, layer_params, perturbation_file, processes, verbose
            )
    else:
        create_constant_vm_file(
            perturbation_file, vm_params["nx"] * vm_params["ny"] * vm_params["nz"]
        )

    if args.fault_damage_zone:
        apply_fault_damage_zone(
            srf_location=args.srf_location,
            vm_params=vm_params,
            pert_f_location=perturbation_file,
            depth_km=args.depth_km,
            max_depth_km=args.max_depth_km,
            width_km=args.width_km,
            max_width_km=args.max_width_km,
            min_damage_velocity=args.max_velocity_drop,
            n_processes=processes,
        )


if __name__ == "__main__":
    main()

import argparse
from os.path import abspath

from qcore.utils import load_yaml

from srf_generation.velocity_model_generation.fault_damage_zone import (
    create_empty_perturbation_file,
    apply_fault_damage_zone,
    add_fault_damage_zone_properties,
)
from srf_generation.velocity_model_generation.generate_3d_velocity_model_perturbation import (
    generate_velocity_model_perturbation_file,
    load_parameter_file,
)


def load_args():
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="A master script to handle all processes"
    )
    parser.add_argument("perturbation_file", type=abspath)
    parser.add_argument("vm_params_location", type=abspath)
    parser.add_argument("-n", "--n_processes", default=1, type=int)

    parser.add_argument("--perturbation", action="store_true")
    parser.add_argument("--fault_damage_zone", action="store_true")

    perturbation_group = parser.add_argument_group()
    perturbation_group.add_argument("--parameter_file", type=abspath)

    fault_damage_zone_group = parser.add_argument_group()
    fault_damage_zone_group.add_argument(
        "--srf_location",
        type=abspath,
        help="Path to the srf. Required if --perturbation is passed, ignored otherwise",
    )
    add_fault_damage_zone_properties(fault_damage_zone_group)

    args = parser.parse_args()

    errors = []

    if args.perturbation:
        if args.parameter_file is None:
            errors.append("If perturbation is given, parameter_file must also be given")

    if args.fault_damage_zone:
        if args.srf_location is None:
            errors.append(
                "If --fault_damage_zone is given, srf_location must also be given"
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

    if args.perturbation:
        common_params, layer_params = load_parameter_file(args.parameter_file)
        generate_velocity_model_perturbation_file(
            common_params, layer_params, perturbation_file, processes
        )
    else:
        create_empty_perturbation_file(
            pert_f_location=perturbation_file, vm_params=vm_params
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

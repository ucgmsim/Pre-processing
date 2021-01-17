import argparse
from multiprocessing.pool import Pool
from os.path import abspath, join, exists
import pandas as pd

from qcore.formats import load_fault_selection_file
from qcore.utils import load_yaml
from qcore.simulation_structure import (
    get_realisation_name,
    get_srf_path,
    get_fault_VM_dir,
)

from srf_generation.velocity_model_generation.fault_damage_zone import (
    add_fault_damage_zone_properties,
    create_empty_perturbation_file,
    apply_fault_damage_zone,
)
from srf_generation.velocity_model_generation.generate_3d_velocity_model_perturbation import (
    generate_velocity_model_perturbation_file_from_config,
    load_parameter_file,
    generate_velocity_model_perturbation_file_from_model,
)


def generate_vm_perturbation(
    args,
    cs_root,
    realisation,
    fault_damage_zone,
    perturbation,
    perturbation_file,
    vm_params,
):
    if exists(perturbation_file):
        continue
    if perturbation:
        if args.model:
            perturbation_model = pd.read_csv(args.model)
            generate_velocity_model_perturbation_file_from_model(
                vm_params, perturbation_model, perturbation_file, 1
            )
        elif args.parameter_file:
            common_params, layer_params = load_parameter_file(args.parameter_file)
            generate_velocity_model_perturbation_file_from_config(
                common_params, layer_params, perturbation_file, 1
            )
    else:
        create_empty_perturbation_file(perturbation_file, vm_params)
    if fault_damage_zone:
        srf_location = get_srf_path(cs_root, realisation)
        apply_fault_damage_zone(
            srf_location,
            vm_params,
            perturbation_file,
            args.depth_km,
            args.max_depth_km,
            args.width_km,
            args.max_width_km,
            args.max_velocity_drop,
            1,
        )


def load_args():
    parser = argparse.ArgumentParser(
        allow_abbrev=False, description="A master script to handle all processes"
    )
    parser.add_argument("cs_root", type=abspath)
    parser.add_argument("fault_selection_file", type=abspath)
    parser.add_argument("-n", "--n_processes", default=1, type=int)

    parser.add_argument("--perturbation", action="store_true")
    parser.add_argument("--fault_damage_zone", action="store_true")

    perturbation_group = parser.add_argument_group()
    perturbation_group.add_argument("--parameter_file", type=abspath)
    perturbation_group.add_argument("--model", type=abspath)

    fault_damage_zone_group = parser.add_argument_group()
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
    args = load_args()

    cs_root = args.cs_root
    faults = load_fault_selection_file(args.fault_selection_file)

    vm_params = load_yaml(args.vm_params_location)
    processes = args.n_processes
    perturbation = args.perturbation
    fault_damage_zone = args.fault_damage_zone

    pert_file_params = []

    for fault, count in faults.items():
        for i in range(1, count + 1):
            realisation = get_realisation_name(fault, i)
            perturbation_file = join(
                get_fault_VM_dir(cs_root, realisation), f"{realisation}.pertb"
            )
            pert_file_params.append(
                (
                    args,
                    cs_root,
                    realisation,
                    fault_damage_zone,
                    perturbation,
                    perturbation_file,
                    vm_params,
                )
            )
    pool = Pool(processes)
    pool.starmap(generate_vm_perturbation, pert_file_params)


if __name__ == "__main__":
    main()

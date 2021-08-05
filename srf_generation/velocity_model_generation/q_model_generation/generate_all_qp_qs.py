"""
To be run manually to generate all Q files for an entire simulation directory at once
"""
import argparse
from multiprocessing import pool
from os.path import abspath

from qcore import simulation_structure
from qcore.formats import load_fault_selection_file
from qcore.utils import load_yaml

from srf_generation.velocity_model_generation.q_model_generation.generate_qp_qs import (
    QP_MODEL_PATH,
    QS_MODEL_PATH,
    MIN_QP,
    MAX_QP,
    MAX_QS,
    MIN_QS,
    MEMORY,
    generate_q_file,
    generate_xy_locations,
)


def generate_qp_qs(args: argparse.Namespace, fault: str):
    vm_params_path = simulation_structure.get_vm_params_path(
        args.cybershake_root, fault
    )
    vm_params = load_yaml(vm_params_path)

    xys = generate_xy_locations(vm_params)

    generate_q_file(
        args.qp_model,
        vm_params["nx"],
        vm_params["ny"],
        vm_params["nz"],
        vm_params["hh"],
        xys,
        simulation_structure.get_fault_qp_file(args.cybershake_root, fault),
        args.min_qp,
        args.max_qp,
        args.useable_ram,
    )
    generate_q_file(
        args.qs_model,
        vm_params["nx"],
        vm_params["ny"],
        vm_params["nz"],
        vm_params["hh"],
        xys,
        simulation_structure.get_fault_qs_file(args.cybershake_root, fault),
        args.min_qs,
        args.max_qs,
        args.useable_ram,
    )


def load_args():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument(
        "cybershake_root", type=abspath, help="Path to the simulation root directory."
    )
    parser.add_argument("fsf", type=abspath, help="Path to the fault selection file.")

    parser.add_argument(
        "-n",
        "--n_processes",
        type=int,
        help="Number of Qp and Qs files to generate at once. At most one process per fault xor event.",
        default=1,
    )

    parser.add_argument(
        "--qp_model",
        type=abspath,
        help="Path to the Qp model to use. Must be in standard VM file binary format. "
        "Must have the same dimensions as the current Eberhart-Phillips model.",
        default=QP_MODEL_PATH,
    )
    parser.add_argument(
        "--qs_model",
        type=abspath,
        help="Path to the Qs model to use. Must be in standard VM file binary format. "
        "Must have the same dimensions as the current Eberhart-Phillips model.",
        default=QS_MODEL_PATH,
    )

    parser.add_argument(
        "--min_qp", type=float, default=MIN_QP, help="Minimum qp value to use."
    )
    parser.add_argument(
        "--max_qp", type=float, default=MAX_QP, help="Maximum qp value to use."
    )
    parser.add_argument(
        "--min_qs", type=float, default=MIN_QS, help="Minimum qs value to use."
    )
    parser.add_argument(
        "--max_qs", type=float, default=MAX_QS, help="Maximum qs value to use."
    )

    parser.add_argument(
        "--useable_ram",
        type=float,
        default=MEMORY,
        help="Maximum available ram to use. In Gb per process.",
    )

    args = parser.parse_args()
    return args


def main():

    args = load_args()
    faults = load_fault_selection_file(args.fsf)

    with pool.Pool(args.n_processes) as worker_pool:
        worker_pool.starmap(generate_qp_qs, [(args, fault) for fault in faults.keys()])


if __name__ == "__main__":
    main()

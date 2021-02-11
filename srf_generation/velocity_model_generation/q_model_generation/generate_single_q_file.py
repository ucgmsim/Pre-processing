"""
To generate a single Q file
"""

import argparse
from os.path import abspath

from qcore.utils import load_yaml

from srf_generation.velocity_model_generation.q_model_generation.generate_qp_qs import (
    generate_xys,
    generate_q_file,
    MAX_QP,
    MIN_QP,
)


def load_args():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument(
        "vm_params_path",
        type=abspath,
        help="Path to the vm_params file to load parameters from.",
    )
    parser.add_argument(
        "model_path",
        type=abspath,
        help="Path to the model to use. Must be in standard VM file binary format. Must have the same dimensions as the current Eberhart-Phillips model.",
    )
    parser.add_argument(
        "outfile_path",
        type=abspath,
        help="Path to the output Q files without suffix.",
    )

    parser.add_argument(
        "--min_q", type=float, default=MIN_QP, help="Minimum qp value to use."
    )
    parser.add_argument(
        "--max_q", type=float, default=MAX_QP, help="Maximum qp value to use."
    )

    args = parser.parse_args()
    return args


def main():
    args = load_args()

    vm_params = load_yaml(args.vm_params_path)

    xys = generate_xys(vm_params)

    generate_q_file(
        args.model_path,
        vm_params["nx"],
        vm_params["ny"],
        vm_params["nz"],
        vm_params["hh"],
        xys,
        args.outfile_path,
        args.min_q,
        args.max_q,
    )


if __name__ == "__main__":
    main()

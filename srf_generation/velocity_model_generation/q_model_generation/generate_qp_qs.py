"""
To generate the Q files for a single realisation
"""

import argparse
import numpy as np
from os.path import abspath
from scipy.interpolate import RegularGridInterpolator

from qcore.geo import gp2ll_multi, ll2xy, gen_mat
from qcore.utils import load_yaml
from qcore.vm_file import VelocityModelFile


Q_LLmat, Q_LLinv = gen_mat(140, 172.9037, -41.7638)
QP_MODEL_PATH = "qpmodel.bin"
QS_MODEL_PATH = "qsmodel.bin"

MAX_QP = np.inf
MIN_QP = 50
MAX_QS = np.inf
MIN_QS = 25

# 2Gb, as is standard on Maui. Memory use is only an estimation, and may not be linear in the number of grid points
MEMORY = 2 * (10 ** 3) ** 3

# fmt: off
XS = [
    -1200, -540, -400, -300, -200, -151, -136, -121, -106,  -96,  -86,  -76,  -68,  -61,  -53,  -46,  -38,  -31,
      -23,  -16,   -8,   -1,    7,   14,   22,   29,   37,   41,   44,   48,   51,   55,   58,   62,   65,   69,
       72,   76,   79,   83,   86,   90,   93,   97,  100,  104,  107,  111,  114,  118,  121,  127,  133,  139,
      145,  152,  161,  178,  205,  230,  260,  300,  350, 1200,
]
YS = [
    -1200, -703, -653, -632, -610, -590, -573, -554, -536, -515, -494, -468, -439, -409, -379, -349, -329, -310,
     -291, -272, -253, -228, -198, -173, -148, -123, -105,  -90,  -80,  -70,  -60,  -50,  -40,  -30,  -20,  -10,
        0,   10,   20,   30,   40,   55,   70,   85,  100,  115,  130,  145,  155,  163,  171,  179,  187,  197,
      224,  269,  313,  358,  388,  418,  448,  477,  507,  537,  567,  600,  639,  686,  736, 1200,
]
ZS = [
      -15,   -1,    1,    3,    5,    8,   15,   23,   30,   34,   38,   42,   48,   55,   65,   85,  105,  130,
      155,  185,  225,  275,  370,  620,  750,
]
# fmt: on


def generate_xys(vm_params: "Dict"):
    # Generate all pairs of xy grid points in the VM domain
    points = [[x, y] for x in range(vm_params["nx"]) for y in range(vm_params["ny"])]
    # Convert these points to lat-lon pairs
    lls = gp2ll_multi(
        points,
        vm_params["MODEL_LAT"],
        vm_params["MODEL_LON"],
        vm_params["MODEL_ROT"],
        vm_params["nx"],
        vm_params["ny"],
        vm_params["hh"],
    )
    # Convert the lat lon pairs for each point into grid points in the model domain
    xys = ll2xy(np.asarray(lls), Q_LLinv)
    return xys


def get_interp(infile: str):
    data = VelocityModelFile(len(XS), len(YS), len(ZS), infile)
    with data as vals:
        q_interp = RegularGridInterpolator(
            (XS, YS, ZS), vals, bounds_error=False, fill_value=None
        )
    return q_interp


def generate_q_file(
    model: str,
    nx: int,
    ny: int,
    nz: int,
    hh: float,
    xys: "List",
    output: str,
    min_val: float = -np.inf,
    max_val: float = np.inf,
):
    # Get the interpolator from the model
    interpolator = get_interp(model)

    # Set the step size from a rudimentary ram use calculation
    step_size = int(np.floor(MEMORY / (240 * nx * ny)))
    if step_size < 1:
        # If we don't have enough ram to do this,
        # then raise an exception and allow the user to increase ram allocation or reduce domain size
        raise MemoryError(
            "Not enough ram to perform this operation. Either allocate more ram, or reduce domain size"
        )

    # Create a new Velocity Model object with the correct dimensions
    qp_data = VelocityModelFile(nx, ny, nz, output)
    qp_data.new()
    # Context manager returns direct access to the underlying array
    with qp_data as data_box:
        # Slice the domain into horizontal layers to reduce memory usage
        for offset in range(0, nz, step_size):
            # For the last layer don't generate more values than required
            if nz - offset < step_size:
                step_size = nz - offset

            # Generate each point in this slice
            locations = [
                [x, y, z * hh]
                for x, y in xys
                for z in range(offset, offset + step_size)
            ]

            # Run the interpolator and
            q_vals = interpolator(locations).reshape((nx, ny, step_size))

            # Ensure we are within reasonable limits
            q_vals = np.minimum(max_val, np.maximum(min_val, q_vals))

            # qp_data.set_values(q_vals)
            data_box[:, :, offset : offset + step_size] = q_vals

        # Save the model before the context manager closes it
        qp_data.save()


def load_args():
    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument(
        "vm_params_path",
        type=abspath,
        help="Path to the vm_params file to load parameters from.",
    )
    parser.add_argument(
        "outfile_prefix",
        type=abspath,
        help="Path to the output Q files without suffix.",
    )

    parser.add_argument(
        "--qp_model",
        type=abspath,
        help="Path to the Qp model to use. Must be in standard VM file binary format. Must have the same dimensions as the current Eberhart-Phillips model.",
        default=QP_MODEL_PATH,
    )
    parser.add_argument(
        "--qs_model",
        type=abspath,
        help="Path to the Qs model to use. Must be in standard VM file binary format. Must have the same dimensions as the current Eberhart-Phillips model.",
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

    args = parser.parse_args()
    return args


def main():
    args = load_args()

    vm_params = load_yaml(args.vm_params_path)

    xys = generate_xys(vm_params)

    print("Generating first q file")
    generate_q_file(
        args.qp_model,
        vm_params["nx"],
        vm_params["ny"],
        vm_params["nz"],
        vm_params["hh"],
        xys,
        f"{args.outfile_prefix}.qp",
        args.min_qp,
        args.max_qp,
    )
    print("Generating second q file")
    generate_q_file(
        args.qs_model,
        vm_params["nx"],
        vm_params["ny"],
        vm_params["nz"],
        vm_params["hh"],
        xys,
        f"{args.outfile_prefix}.qs",
        args.min_qs,
        args.max_qs,
    )
    print("Done")


if __name__ == "__main__":
    main()

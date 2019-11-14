import argparse
from multiprocessing.pool import Pool
from os import remove
from os.path import abspath
import pandas as pd
from subprocess import call, DEVNULL
from tempfile import NamedTemporaryFile

from qcore.binary_version import get_unversioned_bin


FRACTAL3D_BIN = get_unversioned_bin("fractal3D")
LAYER_COMBINE_BIN = get_unversioned_bin("combine_pertb")

PERTURBATION_COMMON_COLUMNS = ["nx", "ny", "grid_spacing"]
PERTURBATION_LAYER_COLUMNS = ["nz", "h_corr", "v_corr" "sigma", "seed"]


def create_perturbated_layer(
    index, nx, ny, nz, grid_spacing, h_corr, v_corr, sigma, seed
):
    layer_file = NamedTemporaryFile(delete=False)
    layer_file.close()
    command = [
        FRACTAL3D_BIN,
        f"nx={nx}",
        f"ny={ny}",
        f"nz={nz}",
        f"h={grid_spacing}",
        f"x_corlen={h_corr}",
        f"y_corlen={h_corr}",
        f"z_corlen={v_corr}",
        f"sigma={sigma}",
        f"seed={seed}",
        f"perturbfile={layer_file.name}",
        "usd_fftw=1",
    ]
    call(command, stderr=DEVNULL)
    return index, nz, layer_file.name


def combine_layers(layers, nx, ny, out_file):
    filelist = NamedTemporaryFile(mode="w", delete=False)
    depth = 0
    with filelist as fl:
        for index, nz, file_name in layers:
            fl.write(f"{file_name} {depth} {depth+nz}\n")
            depth += nz
    filelist.close()
    command = [
        LAYER_COMBINE_BIN,
        f"nx={nx}",
        f"ny={ny}",
        f"nz={nz}",
        f"filelist={filelist.name}",
        f"outfile={out_file}",
    ]
    ret_code = call(command, stderr=DEVNULL)
    if ret_code == 0:
        remove(filelist.name)


def load_parameter_file(parameter_file):
    common_info = pd.read_csv(parameter_file, nrows=1)
    layer_info = pd.read_csv(parameter_file, skiprows=2)
    return common_info.to_dict("records")[0], layer_info.to_dict("index")


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_file", type=abspath)
    parser.add_argument("output_file", type=abspath)
    parser.add_argument("-n", "--n_processes", default=1, type=int)
    args = parser.parse_args()
    return args


def kwarg_map(func, kwargs):
    return func(**kwargs)


def main():
    args = load_args()
    common_params, layer_params = load_parameter_file(args.parameter_file)
    worker_pool = Pool(args.n_processes)
    complete_layer_parameters = [
        (create_perturbated_layer, dict({"index": index}, **l_p, **common_params))
        for index, l_p in layer_params.items()
    ]
    layer_info = sorted(worker_pool.starmap(kwarg_map, complete_layer_parameters))
    combine_layers(
        layer_info, common_params["nx"], common_params["ny"], args.output_file
    )
    for _, _, file in layer_info:
        remove(file)


if __name__ == "__main__":
    main()

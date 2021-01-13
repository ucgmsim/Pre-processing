import argparse
from multiprocessing.pool import Pool
from os import remove
from os.path import abspath
from random import randint

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
        f"nz={depth}",
        f"filelist={filelist.name}",
        f"outfile={out_file}",
    ]
    ret_code = call(command, stderr=DEVNULL)
    if ret_code == 0:
        remove(filelist.name)
    else:
        print(
            f"Non zero return code returned by layer combiner. Error code: {ret_code}"
        )


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
    out_file = args.output_file

    generate_velocity_model_perturbation_file_from_config(
        common_params, layer_params, out_file, args.n_processes
    )


def generate_velocity_model_perturbation_file_from_config(
    common_params, layer_params, out_file, n_processes=1
):
    complete_layer_parameters = [
        (create_perturbated_layer, dict({"index": index}, **l_p, **common_params))
        for index, l_p in layer_params.items()
    ]
    layer_info = sorted(Pool(n_processes).starmap(kwarg_map, complete_layer_parameters))
    combine_layers(layer_info, common_params["nx"], common_params["ny"], out_file)
    for _, _, file in layer_info:
        remove(file)


def generate_velocity_model_perturbation_file_from_model(
    vm_params, perturbation_model, out_file, n_processes=1
):
    max_depth = vm_params["nz"] * vm_params["hh"]
    current_depth = 0

    i = 0

    complete_layer_parameters = []
    while current_depth < max_depth:
        layer = perturbation_model.iloc[i]
        layer_depth = layer["depth"]
        if layer["depth"] + current_depth > max_depth:
            layer_depth = max_depth - current_depth
        nz = int(layer_depth / vm_params["hh"])
        seed = randint(0, 2 ** 31 - 1)

        complete_layer_parameters.append(
            (
                i,
                vm_params["nx"],
                vm_params["ny"],
                nz,
                vm_params["hh"],
                layer["h_corr"],
                layer["v_corr"],
                layer["sigma"],
                seed,
            )
        )

        current_depth += layer_depth
        i += 1

    # Save config files for future use
    pd.DataFrame.from_dict(
        {
            "nx": [vm_params["nx"]],
            "ny": [vm_params["ny"]],
            "grid_spacing": [vm_params["hh"]],
        }
    ).to_csv(f"{out_file}.csv", index=False, mode="w")
    pd.DataFrame(
        complete_layer_parameters,
        columns=["index", "nx", "ny", "hh", "nz", "h_corr", "v_corr", "sigma", "seed"],
    )[["nz", "h_corr", "z_corr", "sigma", "seed"]].to_csv(
        f"{out_file}.csv", index=False, mode="a"
    )

    layer_info = sorted(
        Pool(n_processes).starmap(create_perturbated_layer, complete_layer_parameters)
    )
    combine_layers(layer_info, vm_params["nx"], vm_params["nz"], out_file)
    for _, _, file in layer_info:
        remove(file)


if __name__ == "__main__":
    main()

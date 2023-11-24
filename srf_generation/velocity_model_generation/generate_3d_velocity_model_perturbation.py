import argparse
from tempfile import TemporaryDirectory
from multiprocessing.pool import Pool
from pathlib import Path
from os import remove, makedirs
from os.path import abspath, join
from random import randint
from shutil import rmtree

import pandas as pd
import subprocess

from qcore.binary_version import get_unversioned_bin

DEBUG_OUTPUT_LOCATION = subprocess.DEVNULL

FRACTAL3D_BIN = get_unversioned_bin("fractal3D")
LAYER_COMBINE_BIN = get_unversioned_bin("combine_pertb")

PERTURBATION_COMMON_COLUMNS = ["nx", "ny", "grid_spacing"]
PERTURBATION_LAYER_COLUMNS = ["nz", "h_corr", "v_corr" "sigma", "seed"]

LAYER_GENERATION_ERROR_TEXT = (
    "The process failed at the layer generation, try with n_proc=1"
)


def create_perturbed_layer(
    index,
    nx,
    ny,
    nz,
    grid_spacing,
    h_corr,
    v_corr,
    sigma,
    seed,
    temp_dir=Path(abspath(".")),
):
    layer_file = temp_dir / f"layer{index}.pertb"
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
        f"perturbfile={str(layer_file)}",
        "use_fftw=1",
    ]
    print(command)
    subprocess.call(command, stderr=DEBUG_OUTPUT_LOCATION)
    return index, nz, layer_file


def combine_layers(layers, nx, ny, out_file, temp_dir=Path(".").resolve()):
    filelist = temp_dir / "filelist.txt"
    depth = 0
    with open(filelist, "w") as fl:
        for index, nz, file_name in layers:
            fl.write(f"{file_name} {depth} {depth+nz}\n")
            depth += nz
    command = [
        LAYER_COMBINE_BIN,
        f"nx={nx}",
        f"ny={ny}",
        f"nz={depth}",
        f"filelist={filelist}",
        f"outfile={out_file}",
    ]
    print(command)
    ret_code = subprocess.call(command, stderr=DEBUG_OUTPUT_LOCATION)
    if ret_code == 0:
        remove(filelist)
    else:
        print(
            f"Non zero return code returned by layer combiner. Error code: {ret_code}"
        )
        raise AssertionError(
            "ERROR Layer Combiner error code was non-zero. (See stdout for more information)"
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
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Output debug info to std_err"
    )
    args = parser.parse_args()
    return args


def kwarg_map(func, kwargs):
    return func(**kwargs)


def main():
    args = load_args()

    common_params, layer_params = load_parameter_file(args.parameter_file)
    out_file = args.output_file

    generate_velocity_model_perturbation_file_from_config(
        common_params, layer_params, out_file, args.n_processes, args.verbose
    )


def generate_velocity_model_perturbation_file_from_config(
    common_params, layer_params, out_file, n_processes=1, verbose=False
):
    set_verbose(verbose)

    temp_dir = TemporaryDirectory(dir=Path(out_file).parent)
    temp_dir_path = Path(temp_dir.name).resolve()

    complete_layer_parameters = [
        (
            create_perturbed_layer,
            dict({"index": index, "temp_dir": temp_dir_path}, **l_p, **common_params),
        )
        for index, l_p in layer_params.items()
    ]

    if n_processes == 1:
        layer_info = sorted([kwarg_map(*layer) for layer in complete_layer_parameters])
    else:
        try:
            with Pool(n_processes) as pool:
                layer_info = sorted(pool.starmap(kwarg_map, complete_layer_parameters))
        except AssertionError:
            print("The process failed at the layer generation, try with n_proc=1")
            raise
    combine_layers(layer_info, common_params["nx"], common_params["ny"], out_file, temp_dir_path)
    for _, _, file in layer_info:
        remove(file)
    temp_dir.cleanup()


def set_verbose(verbose):
    if verbose:
        global DEBUG_OUTPUT_LOCATION
        DEBUG_OUTPUT_LOCATION = None


def generate_velocity_model_perturbation_file_from_model(
    vm_params, perturbation_model, out_file, n_processes=1, verbose=False
):
    set_verbose(verbose)
    max_depth = vm_params["nz"] * vm_params["hh"]
    current_depth = 0

    i = 0

    # Uses a temp dir in the current folder structure due to memory issues if using an actual temp directory
    with TemporaryDirectory(dir=Path(out_file).parent) as temp_dir:
        temp_dir_path = Path(temp_dir).resolve()
        complete_layer_parameters = []
        while current_depth < max_depth:
            layer = perturbation_model.iloc[i]
            layer_depth = layer["depth"]
            if layer["depth"] + current_depth > max_depth:
                layer_depth = max_depth - current_depth
            nz = int(layer_depth / vm_params["hh"])
            seed = randint(0, 2**31 - 1)

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
                    temp_dir_path,
                )
            )
            # Adds the actual depth of the new layer to the current depth.
            # In some cases this may not be the whole model layer due to rounding.
            current_depth += nz * vm_params["hh"]
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
            columns=[
                "index",
                "nx",
                "ny",
                "nz",
                "hh",
                "h_corr",
                "v_corr",
                "sigma",
                "seed",
                "temp_dir",
            ],
        )[["nz", "h_corr", "v_corr", "sigma", "seed"]].to_csv(
            f"{out_file}.csv", index=False, mode="a"
        )
        if n_processes == 1:
            layer_info = sorted(
                [create_perturbed_layer(*layer) for layer in complete_layer_parameters]
            )
        else:
            try:
                with Pool(n_processes) as pool:
                    layer_info = sorted(
                        pool.starmap(create_perturbed_layer, complete_layer_parameters)
                    )
            except AssertionError:
                print(LAYER_GENERATION_ERROR_TEXT)
                raise
        combine_layers(
            layer_info, vm_params["nx"], vm_params["ny"], out_file, temp_dir_path
        )


if __name__ == "__main__":
    main()

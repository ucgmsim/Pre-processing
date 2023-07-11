import argparse
from os.path import abspath
from random import randint

import pandas as pd

from qcore.utils import load_yaml


def generate_pertb_config(vm_params, perturbation_model: pd.DataFrame):
    max_depth = vm_params["nz"] * vm_params["hh"]
    bottom_depth = 0

    i = 0

    layers = []

    while bottom_depth < max_depth:
        layer = perturbation_model.iloc[i]
        layer_depth = layer["depth"]
        if layer["depth"] + bottom_depth > max_depth:
            layer_depth = max_depth - bottom_depth
        nz = int(layer_depth / vm_params["hh"])
        seed = randint(0, 2**31 - 1)
        layers.append(
            {
                "nz": nz,
                "h_corr": layer["h_corr"],
                "v_corr": layer["v_corr"],
                "sigma": layer["sigma"],
                "seed": seed,
            }
        )

        bottom_depth += layer_depth
        i += 1

    layer_config = pd.DataFrame(layers)

    return generate_config_header(vm_params), layer_config


def generate_config_header(vm_params):
    config_header = {
        "nx": [vm_params["nx"]],
        "ny": [vm_params["ny"]],
        "grid_spacing": [vm_params["hh"]],
    }
    return pd.DataFrame.from_dict(config_header)


def write_pertb_config(
    config_header: pd.DataFrame, config_layers: pd.DataFrame, out_file
):
    config_header.to_csv(out_file, index=False, mode="w")
    config_layers.to_csv(out_file, index=False, mode="a")


def load_args():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("vm_params", type=abspath)
    parser.add_argument("perturbation_model", type=abspath)
    parser.add_argument("out_file", type=abspath)
    args = parser.parse_args()
    return args


def main():
    args = load_args()
    vm_params = load_yaml(args.vm_params)
    perturbation_model = pd.read_csv(args.perturbation_model)
    config_header, config_layers = generate_pertb_config(vm_params, perturbation_model)
    write_pertb_config(config_header, config_layers, args.out_file)


if __name__ == "__main__":
    main()

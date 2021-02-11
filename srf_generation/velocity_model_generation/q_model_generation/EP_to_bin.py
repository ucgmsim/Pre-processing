"""Script to generate a Qp/Qs binary model file from an Eberhart-Phillips style ascii model"""
import argparse
import numpy as np
from os.path import abspath
import pandas as pd

from qcore.vm_file import VelocityModelFile


def load_args():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("infile", type=abspath)
    parser.add_argument("outfile", type=abspath)
    args = parser.parse_args()
    return args


def main():
    args = load_args()

    # load model data
    data = pd.read_csv(args.infile, sep="\s+", skiprows=1, header=0, index_col=None)
    data.columns = ["vals", "SF", "x", "y", "z", "lat", "lon"]
    data.sort_values(["x", "y", "z"], axis=0, inplace=True)
    data.reset_index(inplace=True, drop=True)
    xs = data["x"].unique()
    ys = data["y"].unique()
    zs = data["z"].unique()
    xs_count = len(xs)
    ys_count = len(ys)
    zs_count = len(zs)
    vals = data["vals"].values.reshape((xs_count, ys_count, zs_count))

    # Save model data
    q_model = VelocityModelFile(xs_count, ys_count, zs_count, args.outfile)
    q_model.new()
    with q_model:
        q_model.set_values(vals)
        q_model.save()

    with open(f"{args.outfile}.meta", "w") as metaf:
        metaf.write(f"Xs: {xs}\n")
        metaf.write(f"Ys: {ys}\n")
        metaf.write(f"Zs: {zs}\n")
    # Present user with model metadata
    print(f"File saved to {args.outfile}")
    print(f"Metadata saved to {args.outfile}.meta")
    print(f"x axis values: {xs}")
    print(f"y axis values: {ys}")
    print(f"z axis values: {zs}")


if __name__ == "__main__":
    main()

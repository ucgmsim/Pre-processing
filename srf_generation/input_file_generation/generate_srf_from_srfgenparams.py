import argparse
import glob
from multiprocessing.pool import Pool
from os import path

import pandas as pd

from qcore import simulation_structure

from srf_generation.input_file_generation.srfgenparams_to_srf import create_ps_srf, generate_sim_params_yaml


def process_srfgenparams_file(cybershake_root, srfgenparams_file):
    rel_df: pd.DataFrame = pd.read_csv(srfgenparams_file)
    realisation = rel_df.to_dict(orient="records")[0]
    if realisation["type"] == 1:
        create_ps_srf(cybershake_root, realisation)
    generate_sim_params_yaml(cybershake_root, realisation)


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("cybershake_root", type=path.abspath)
    parser.add_argument("-n", "--n_processes", default=1, type=int)
    return parser.parse_args()


def main():
    args = load_args()

    srfgenparams_path = simulation_structure.get_srf_path(
        args.cybershake_root, "*"
    ).replace(".srf", ".csv")
    srfgenparams_files = glob.glob(srfgenparams_path)

    worker_pool = Pool(args.n_processes)

    worker_pool.starmap(
        process_srfgenparams_file,
        [(args.cybershake_root, filename) for filename in srfgenparams_files],
    )


if __name__ == "__main__":
    main()

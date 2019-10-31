import argparse
import glob
from multiprocessing.pool import Pool
from os import path

import pandas as pd

from qcore import simulation_structure

from srf_generation.input_file_generation.srfgenparams_to_srf import (
    create_ps_srf,
    generate_sim_params_yaml,
    create_info_file,
)


def process_srfgenparams_file(cybershake_root, srfgenparams_file):
    rel_df: pd.DataFrame = pd.read_csv(srfgenparams_file)
    realisation = rel_df.to_dict(orient="records")[0]
    if realisation["type"] == 1:
        create_ps_srf(cybershake_root, realisation)
    sim_params_file = simulation_structure.get_source_params_path(
        cybershake_root, realisation["name"]
    )
    generate_sim_params_yaml(sim_params_file, realisation)


def process_common_srfgenparams_file(cybershake_root, srfgenparams_file):
    rel_df: pd.DataFrame = pd.read_csv(srfgenparams_file)
    realisation = rel_df.to_dict(orient="records")[0]
    srf_file = simulation_structure.get_srf_path(
        cybershake_root,
        simulation_structure.get_realisation_name(realisation["name"], 1),
    )
    info_filename = srfgenparams_file.replace(".csv", ".info")
    if realisation["type"] == 1:
        create_info_file(
            srf_file=srf_file,
            srf_type=realisation["type"],
            mag=realisation["magnitude"],
            rake=realisation["rake"],
            dt=realisation.get("dt", 0.005),
            vs=realisation.get("vs", 3.2),
            rho=realisation.get("rho", 2.44),
            file_name=info_filename,
            lon=realisation["longitude"],
            lat=realisation["latitude"],
            centroid_depth=realisation["depth"],
        )
    else:
        raise ValueError(
            f"Type {realisation['type']} faults are not currently supported. "
            f"Contact the software team if you believe this is an error."
        )

    sim_params_file = simulation_structure.get_source_params_path(
        cybershake_root, realisation["name"]
    )
    generate_sim_params_yaml(sim_params_file, realisation)


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cybershake_root", type=path.abspath, default=path.abspath("."))
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

    srfgenparams_path = path.join(
        simulation_structure.get_sources_dir(args.cybershake_root), "*", "*.csv"
    )
    srfgenparams_files = glob.glob(srfgenparams_path)

    worker_pool = Pool(args.n_processes)

    worker_pool.starmap(
        process_common_srfgenparams_file,
        [(args.cybershake_root, filename) for filename in srfgenparams_files],
    )


if __name__ == "__main__":
    main()

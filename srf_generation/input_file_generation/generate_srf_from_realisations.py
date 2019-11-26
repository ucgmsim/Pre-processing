import argparse
import glob
from logging import Logger
from multiprocessing.pool import Pool
from os import path
from typing import Dict

import pandas as pd

from qcore import simulation_structure, qclogging
from qcore.qclogging import NOPRINTCRITICAL, add_general_file_handler

from srf_generation.input_file_generation.realisation_to_srf import (
    create_ps_srf,
    generate_sim_params_yaml,
    create_info_file,
    create_ps_ff_srf,
)


def process_realisation_file(
    cybershake_root: str, realisation_file: str, logger_name: str
):

    primary_logger = qclogging.get_logger(logger_name)
    primary_logger.debug(f"Getting realisation values for {realisation_file}")

    rel_df: pd.DataFrame = pd.read_csv(realisation_file, dtype={"name": str})
    realisation = rel_df.to_dict(orient="records")[0]

    realisation_logger = qclogging.get_realisation_logger(
        primary_logger, realisation["name"]
    )
    realisation_logger.debug(f"Got realisation file info: {realisation}")
    stoch_file = simulation_structure.get_stoch_path(
        cybershake_root, realisation["name"]
    )

    if realisation["type"] == 1:
        realisation_logger.debug(
            "Realisation is of type 1, generating srf and related files"
        )
        create_ps_srf(
            realisation_file, realisation, stoch_file, logger=realisation_logger
        )
    elif realisation["type"] == 2:
        realisation_logger.debug(
            "Realisation is of type 2, generating srf and related files"
        )
        create_ps_ff_srf(
            realisation_file, realisation, stoch_file, logger=realisation_logger
        )
    else:
        message = f"The realisation type is not valid: {realisation['type']}, aborting."
        realisation_logger.log(NOPRINTCRITICAL, message)
        raise ValueError(message)

    sim_params_file = simulation_structure.get_source_params_path(
        cybershake_root, realisation["name"]
    )
    generate_sim_params_yaml(sim_params_file, realisation, realisation_logger)


def process_common_realisation_file(
    cybershake_root, realisation_file, logger: Logger = qclogging.get_basic_logger()
):
    rel_df: pd.DataFrame = pd.read_csv(realisation_file, dtype={"name": str})
    realisation: Dict = rel_df.to_dict(orient="records")[0]
    srf_file = simulation_structure.get_srf_path(
        cybershake_root,
        simulation_structure.get_realisation_name(realisation["name"], 1),
    )
    info_filename = realisation_file.replace(".csv", ".info")
    if realisation["type"] == 1:
        create_info_file(
            srf_file=srf_file,
            srf_type=realisation.get("type"),
            mag=realisation.get("magnitude"),
            rake=realisation.get("rake"),
            dt=realisation.get("dt", 0.005),
            vs=realisation.get("vs", 3.2),
            rho=realisation.get("rho", 2.44),
            file_name=info_filename,
            lon=realisation.get("longitude"),
            lat=realisation.get("latitude"),
            centroid_depth=realisation.get("depth"),
        )
    elif realisation["type"] == 2:
        create_info_file(
            srf_file,
            2,
            realisation.get("magnitude"),
            realisation.get("rake"),
            realisation.get("dt", 0.005),
            centroid_depth=realisation.get("depth"),
            lon=None,
            lat=None,
            tect_type=None,
            dip_dir=None,
            mwsr=realisation.get("depth"),
            shypo=realisation.get("shypo") + 0.5 * realisation.get("flen"),
            dhypo=realisation.get("dhypo"),
            vm=realisation.get("vel_mod_1d"),
            logger=logger,
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
    parser.add_argument(
        "--cybershake_root", type=path.abspath, default=path.abspath(".")
    )
    parser.add_argument("-n", "--n_processes", default=1, type=int)
    return parser.parse_args()


def main():
    primary_logger = qclogging.get_logger("generate_srf_from_realisations")
    args = load_args()
    log_file = path.join(args.cybershake_root, primary_logger.name + "_log.txt")
    add_general_file_handler(primary_logger, log_file)
    primary_logger.debug(f"Args parsed: {args}")

    realisations_path = simulation_structure.get_srf_path(
        args.cybershake_root, "*"
    ).replace(".srf", ".csv")
    primary_logger.debug(
        f"Checking for realisation files that match the following path format: {realisations_path}"
    )
    realisation_files = glob.glob(realisations_path)
    primary_logger.debug(f"Got the following realisation files: {realisation_files}")

    worker_pool = Pool(args.n_processes)

    worker_pool.starmap(
        process_realisation_file,
        [
            (args.cybershake_root, filename, primary_logger.name)
            for filename in realisation_files
        ],
    )

    realisations_path = path.join(
        simulation_structure.get_sources_dir(args.cybershake_root), "*", "*.csv"
    )
    realisation_files = glob.glob(realisations_path)

    worker_pool.starmap(
        process_common_realisation_file,
        [(args.cybershake_root, filename) for filename in realisation_files],
    )


if __name__ == "__main__":
    main()

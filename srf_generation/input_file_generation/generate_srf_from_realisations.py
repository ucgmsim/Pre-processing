import argparse
import glob
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
    create_multi_plane_srf,
)
from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)


def process_realisation_file(
    cybershake_root: str, realisation_file: str, checkpointing: bool, logger_name: str
):

    primary_logger = qclogging.get_logger(logger_name)
    primary_logger.debug(f"Getting realisation values for {realisation_file}")

    rel_df: pd.DataFrame = pd.read_csv(realisation_file, dtype={"name": str})
    realisation = rel_df.to_dict(orient="records")[0]

    realisation_logger = qclogging.get_realisation_logger(
        primary_logger, realisation["name"]
    )

    if checkpointing and path.isfile(realisation_file.replace(".csv", ".info")):
        realisation_logger.debug(
            f"Info file for realisation {realisation['name']} already exists, continuing"
        )
        return

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
    elif realisation["type"] == 4:
        realisation_logger.debug(
            "Realisation is of type 4, generating srf and related files"
        )
        create_multi_plane_srf(
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
    cybershake_root, realisation_file, checkpointing: bool, logger_name: str
):
    logger = qclogging.get_logger(logger_name)
    rel_df: pd.DataFrame = pd.read_csv(realisation_file, dtype={"name": str})
    realisation: Dict = rel_df.to_dict(orient="records")[0]
    srf_file = simulation_structure.get_srf_path(
        cybershake_root,
        simulation_structure.get_realisation_name(realisation["name"], 1),
    )
    info_filename = realisation_file.replace(".csv", ".info")
    if checkpointing and path.isfile(info_filename):
        logger.debug(
            f"Common info file for realisation {realisation['name']} already exists, continuing"
        )
        return
    if realisation["type"] == 1:
        create_info_file(
            srf_file=srf_file,
            srf_type=1,
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
            file_name=info_filename,
            centroid_depth=realisation.get("depth"),
            lon=None,
            lat=None,
            tect_type=realisation.get("tect_type", None),
            dip_dir=None,
            mwsr=realisation.get("mwsr"),
            shypo=realisation.get("shypo") + 0.5 * realisation.get("flen"),
            dhypo=realisation.get("dhypo"),
            vm=realisation.get("srf_vel_mod_1d", DEFAULT_1D_VELOCITY_MODEL_PATH),
            logger=logger,
        )
    elif realisation["type"] == 4:
        create_info_file(
            srf_file,
            4,
            realisation.get("magnitude"),
            realisation.get("rake"),
            realisation.get("dt", 0.005),
            file_name=info_filename,
            tect_type=realisation.get("tect_type"),
            dip_dir=realisation.get("dip_dir"),
            shypo=[
                realisation.get("shypo") + 0.5 * realisation.get(f"length_subfault_{i}")
                for i in range(realisation.get("plane_count"))
            ],
            dhypo=realisation.get("dhypo"),
            vm=realisation.get("srf_vel_mod_1d", DEFAULT_1D_VELOCITY_MODEL_PATH),
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
    parser.add_argument(
        "-c",
        "--checkpointing",
        action="store_true",
        help="Activate checkpointing for srf/stoch/info generation. Any realisations with an info file already exists will not be regenerated. Default is to regenerate.",
    )
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
    realisation_files = sorted(glob.glob(realisations_path))
    primary_logger.debug(f"Got the following realisation files: {realisation_files}")

    worker_pool = Pool(args.n_processes)

    worker_pool.starmap(
        process_realisation_file,
        [
            (args.cybershake_root, filename, args.checkpointing, primary_logger.name)
            for filename in realisation_files
        ],
    )

    realisations_path = path.join(
        simulation_structure.get_sources_dir(args.cybershake_root), "*", "*.csv"
    )
    realisation_files = glob.glob(realisations_path)

    worker_pool.starmap(
        process_common_realisation_file,
        [
            (args.cybershake_root, filename, args.checkpointing, primary_logger.name)
            for filename in realisation_files
        ],
    )


if __name__ == "__main__":
    main()

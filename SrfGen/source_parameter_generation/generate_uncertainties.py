import argparse
from multiprocessing import pool
from os import makedirs
from os.path import abspath, isfile, isdir, dirname
import pandas as pd
import yaml

from qcore.formats import load_fault_selection_file
from qcore.simulation_structure import get_realisation_name, get_srf_path

from SrfGen.source_parameter_generation.uncertainties.versions import PERTURBATORS


def load_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("version", type=str)
    parser.add_argument("station_file")
    parser.add_argument("gcmt_file")
    parser.add_argument("-n", "--n_processes", default=1)
    parser.add_argument("-c", "--cybershake_root", type=str, default=".")

    args = parser.parse_args()

    args.station_file = abspath(args.station_file)
    args.gcmt_file = abspath(args.gcmt_file)

    errors = []

    if args.version not in PERTURBATORS.keys():
        errors.append(f"Specified version not found: {args.version}")
    if not isfile(args.station_file):
        errors.append(f"Specified station file not found: {args.station_file}")
    if not isfile(args.gcmt_file):
        errors.append(f"Specified gcmt file not found: {args.gcmt_file}")

    if errors:
        message = "At least one error was detected when verifying arguments:\n"+"\n".join(errors)
        raise ValueError(message)

    makedirs(args.cybershake_root, exist_ok=True)

    return args


def generate_unceratinties(data: pd.Series, realisation_count: int, cybershake_root: str, perturbator_version: str):
    for i in range(realisation_count):
        rel_name = get_realisation_name(data.PublicID, i)
        perturbed_realisation = PERTURBATORS[perturbator_version](data)
        file_name = get_srf_path(cybershake_root, rel_name).replace('.srf', '.yaml')
        makedirs(dirname(file_name))
        with open(file_name) as params_file:
            yaml.dump(perturbed_realisation, params_file)


def main():
    args = load_args()
    faults = load_fault_selection_file(args.station_file)
    gcmt_data = pd.read_csv(
        args.gcmt_file
    )

    gcmt_lines = gcmt_data[gcmt_data['PublicID'].isin(faults.keys())].itertuples(index=False)
    missing_faults = [fault for fault in faults.keys() if fault not in gcmt_lines['PublicID']]
    if missing_faults:
        raise ValueError("Some fault(s) in the fault selection were not found in the gcmt file: {}".format(missing_faults))

    messages = [
        (gcmt_lines[i], faults[gcmt_lines[i].PublicID], args.cybershake_root, args.version)
        for i in range(len(faults))
    ]

    worker_pool = pool.Pool(processes=args.n_processes)
    worker_pool.starmap(generate_unceratinties, messages)


if __name__ == '__main__':
    main()

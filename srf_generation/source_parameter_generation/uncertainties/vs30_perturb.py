
import argparse

import numpy as np
import pandas as pd
from pathlib import Path

from qcore.uncertainties import distributions
from srf_generation.source_parameter_generation import gcmt_to_realisation as gcmt2rel



def perturb_vs30(vs30_data, measurement_uncertainty=0.06, site_model_uncertainty=0.3):
    total_vs30_uncertainty = np.sqrt(
        (vs30_data["sigma"].values) ** 2 + measurement_uncertainty**2
    )

    assert vs30_data is not None
    perturbed_vs30 = vs30_data.copy(deep=True)
    perturbed_vs30["vs30"] = distributions.truncated_log_normal(vs30_data["median"].values, total_vs30_uncertainty, 2)
    return perturbed_vs30
    
    
if __name__ == '__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description='Perturb Vs30 for N realisations')

    # Add an argument for the Vs30 file path
    parser.add_argument("vs30_median", type=Path, help="The path to a file containing median VS30s.")
    parser.add_argument('n_rels', type=int, help='Number of realisations')
    parser.add_argument(
            "outdir", type=Path, help="The path for output n_rels vs30 files"
    )
    parser.add_argument(
            "--vs30_sigma", type=Path, help="The path to a file containing VS30 sigmas."
    )

    # Parse the command-line arguments
    args = parser.parse_args()
    args.vs30_median= args.vs30_median.resolve()
    args.outdir = args.outdir.resolve()

    if args.vs30_sigma:
        args.vs30_sigma = args.vs30_sigma.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)

    vs30_df = gcmt2rel.load_vs30_median_sigma(args.vs30_median, args.vs30_sigma)
    
    vs30_median_name=args.vs30_median.stem
    for i in range(args.n_rels):
        perturbed_vs30 = perturb_vs30(vs30_df)
        perturbed_vs30.to_csv(args.outdir/f"{vs30_median_name}_REL{i+1:02d}.vs30", columns=["vs30"], sep=" ", index=True, header=False)

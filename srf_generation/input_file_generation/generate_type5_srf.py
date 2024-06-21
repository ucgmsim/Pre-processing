#!usr/bin/env python
"""
Type-5 SRF file generation.

This script facilitates the generation of type-5 SRF files from realisation
specifications (YAML files). It generates fault geometries, generates SRF
files from fault geometries, and stitches together these files in the order
of rupture propagation.

Usage
-----
To generate a type-5 SRF file from a realisation specification:

```
$ python srf_generation.py path/to/realisation.yaml output_directory
```
"""

import functools
import multiprocessing
import subprocess
from pathlib import Path
from typing import Annotated

import numpy as np
import pandas as pd
import srf
import typer
import yaml

from qcore import binary_version, coordinates, grid, gsf
from srf_generation import realisation
from srf_generation.realisation import Realisation, RealisationFault
from srf_generation.source_parameter_generation import uncertainties

FAULTSEG2GSFDIPDIR = "fault_seg2gsf_dipdir"
SRF2STOCH = "srf2stoch"


def generate_fault_gsf(
    fault: RealisationFault,
    gsf_output_directory: Path,
    subdivision_resolution: float,
):
    """Write the fault geometry of a fault to a GSF file.

    Parameters
    ----------
    fault : RealisationFault
        The fault to write.
    gsf_output_directory : Path
        The directory to output the GSF file to.
    subdivision_resolution : float
        The geometry resolution.
    """
    gsf_output_filepath = (
        gsf_output_directory / f"{realisation.normalise_name(fault.name)}.gsf"
    )
    gsf_df = pd.DataFrame(
        [
            {
                "strike": plane.strike,
                "dip": plane.dip,
                "length": plane.length,
                "width": plane.width,
                "rake": plane.rake,
                "meshgrid": grid.coordinate_meshgrid(
                    plane.corners[0],
                    plane.corners[1],
                    plane.corners[-1],
                    subdivision_resolution * 1000,
                ),
            }
            for plane in fault.planes
        ]
    )
    gsf.write_fault_to_gsf_file(
        gsf_output_filepath, gsf_df, subdivision_resolution * 1000
    )

    return gsf_output_filepath


def create_stoch(
    srf_file: Path,
    stoch_file: Path,
):
    """Create a stoch file from a SRF file.

    Parameters
    ----------
    srf_file : Path
        The filepath of the SRF file.
    stoch_file : Path
        The filepath to output the stoch file to.
    """

    dx = 2.0
    dy = 2.0
    srf2stoch = binary_version.get_unversioned_bin(SRF2STOCH)
    command = [
        srf2stoch,
        f"dx={dx}",
        f"dy={dy}",
        f"infile={str(srf_file)}",
        f"outfile={str(stoch_file)}",
    ]
    subprocess.run(command, stderr=subprocess.PIPE, check=True)


def generate_fault_srf(
    fault: RealisationFault,
    realisation: Realisation,
    output_directory: Path,
    subdivision_resolution: float,
):
    """Generate an SRF file for a given fault.

    Parameters
    ----------
    fault : RealisationFault
        The fault to generate the SRF file for.
    realisation : Realisation
        The realisation the fault belongs to.
    output_directory : Path
        The output directory.
    subdivision_resolution : float
        The geometry resolution.
    """
    gsf_output_directory = output_directory / "gsf"
    gsf_file_path = generate_fault_gsf(
        fault, gsf_output_directory, subdivision_resolution
    )

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)

    nx = sum(
        grid.gridpoint_count_in_length(plane.length_m, subdivision_resolution * 1000)
        for plane in fault.planes
    )
    ny = grid.gridpoint_count_in_length(
        fault.planes[0].width_m, subdivision_resolution * 1000
    )
    genslip_cmd = [
        genslip_bin,
        "read_erf=0",
        "write_srf=1",
        "read_gsf=1",
        "write_gsf=0",
        f"infile={gsf_file_path}",
        f"mag={fault.magnitude}",
        f"nstk={nx}",
        f"ndip={ny}",
        "ns=1",
        "nh=1",
        f"seed={realisation.genslip_seed}",
        f"velfile={realisation.velocity_model}",
        f"shypo={fault.shyp or 0}",
        f"dhypo={fault.dhyp or 0}",
        f"dt={realisation.dt}",
        "plane_header=1",
        "srf_version=1.0",
        "seg_delay={0}",
        "rvfac_seg=-1",
        "gwid=-1",
        "side_taper=0.02",
        "bot_taper=0.02",
        "top_taper=0.0",
        "rup_delay=0",
        "alpha_rough=0.0",
    ]

    if fault.tect_type == "SUBDUCTION_INTERFACE":
        genslip_cmd.extend(
            [
                "kmodel=-1",
                "xmag_exp=0.5",
                "ymag_exp=0.5",
                "kx_corner=2.5482",
                "ky_corner=2.3882",
                "tsfac_slope=-0.5",
                "tsfac_bzero=-0.1",
                "risetime_coef=1.95",
            ]
        )
    srf_file_path = output_directory / (fault.name + ".srf")
    with open(srf_file_path, "w", encoding="utf-8") as srf_file_handle:
        subprocess.run(
            genslip_cmd, stdout=srf_file_handle, stderr=subprocess.PIPE, check=True
        )


def stitch_srf_files(realisation_obj: Realisation, output_directory: Path) -> Path:
    """Stitch SRF files together in the order of rupture propogation.

    Parameters
    ----------
    realisation_obj : Realisation
        The realisation containing the faults and rupture propogation order.
    output_directory : Path
        The output directory containing fault SRF files.

    Returns
    -------
    Path
        The path to the stitched together SRF file.
    """
    srf_output_filepath = output_directory / f"{realisation_obj.name}.srf"
    with open(srf_output_filepath, "w", encoding="utf-8") as srf_file_output:
        fault_points = {}
        header = []

        srf.write_version(srf_file_output)

        for fault in realisation.topologically_sorted_faults(
            list(realisation_obj.faults.values())
        ):
            with open(
                output_directory / (fault.name + ".srf"), "r", encoding="utf-8"
            ) as fault_srf_file:
                srf.read_version(fault_srf_file)
                fault_header = srf.read_srf_headers(fault_srf_file)
                if fault.parent:
                    for plane in fault_header:
                        # The value of -999, -999 is used in the SRF spec to say
                        # "no hypocentre for this segment".
                        plane.shyp = -999
                        plane.dhyp = -999
                header.extend(fault_header)
                point_count = srf.read_points_count(fault_srf_file)
                fault_points[fault.name] = srf.read_srf_n_points(
                    point_count, fault_srf_file
                )

        srf.write_srf_header(srf_file_output, header)
        srf.write_point_count(
            srf_file_output, sum(len(points) for points in fault_points.values())
        )
        for fault in realisation.topologically_sorted_faults(
            list(realisation_obj.faults.values())
        ):
            t_delay = 0
            if fault.parent:
                # find closest grid point to the jump location
                # compute the time delay as equal to the tinit of this point (for now)
                parent = fault.parent
                parent_coords = coordinates.wgs_depth_to_nztm(
                    np.array(fault.parent_jump_coords)
                )
                parent_fault_points = fault_points[parent.name]
                grid_points = coordinates.wgs_depth_to_nztm(
                    np.array(
                        [
                            [point.lat, point.lon, point.dep * 1000]
                            for point in parent_fault_points
                        ]
                    )
                )
                jump_index = int(
                    np.argmin(np.sum(np.square(grid_points - parent_coords), axis=1))
                )
                t_delay = parent_fault_points[jump_index].tinit
            cur_fault_points = fault_points[fault.name]
            for point in cur_fault_points:
                point.tinit += t_delay
                srf.write_srf_point(srf_file_output, point)

        return srf_output_filepath


def generate_fault_srfs_parallel(
    realisation: Realisation,
    output_directory: Path,
    subdivision_resolution: float,
):
    """Generate fault SRF files in parallel.

    Parameters
    ----------
    realisation : Realisation
        The realisation to generate fault SRF files for.
    output_directory : Path
        The directory to output the fault SRF files.
    subdivision_resolution : float
        The geometry resolution.
    """
    # need to do this before multiprocessing because of race conditions
    gsf_directory = output_directory / "gsf"
    gsf_directory.mkdir(exist_ok=True)

    with multiprocessing.Pool() as worker_pool:
        worker_pool.map(
            functools.partial(
                generate_fault_srf,
                output_directory=output_directory,
                realisation=realisation,
                subdivision_resolution=subdivision_resolution,
            ),
            realisation.faults.values(),
        )


def generate_sim_params_yaml(
    sim_params_filepath: Path,
):
    """Write simulation parameters to a yaml file.

    Parameters
    ----------
    sim_params_file : str
        The filepath to save the simulation parameters.
    """
    # stolen from realisation_to_srf
    sim_params = {}

    with open(sim_params_filepath, "w", encoding="utf-8") as sim_params_file_handle:
        yaml.safe_dump(sim_params, sim_params_file_handle)


def main(
    realisation_filepath: Annotated[
        Path,
        typer.Argument(
            exists=True,
            readable=True,
            help="The filepath of the YAML file containing the realisation data.",
            dir_okay=False,
        ),
    ],
    output_directory: Annotated[
        Path,
        typer.Argument(
            exists=True,
            writable=True,
            file_okay=False,
            help="The output directory path for SRF files.",
        ),
    ],
    subdivision_resolution: Annotated[
        float, typer.Option(help="Geometry resolution (in km)", min=0)
    ] = 0.1,
):
    """Generate a type-5 SRF file from a given realisation specification."""
    realisation_parameters = realisation.read_realisation(realisation_filepath)
    generate_fault_srfs_parallel(
        realisation_parameters, output_directory, subdivision_resolution
    )
    srf_output_filepath = stitch_srf_files(realisation_parameters, output_directory)
    create_stoch(
        srf_output_filepath, output_directory / f"{realisation_parameters.name}.stoch"
    )
    generate_sim_params_yaml(
        output_directory / realisation_filepath.name.removeprefix("type5_")
    )


if __name__ == "__main__":
    typer.run(main)

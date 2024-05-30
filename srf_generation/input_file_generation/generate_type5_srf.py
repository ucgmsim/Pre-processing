import multiprocessing
import re
import shlex
import subprocess
from pathlib import Path
from typing import Annotated

import numpy as np
import pandas as pd
import pyproj
import srf
import typer
from qcore import binary_version, gsf

from srf_generation import realisation
from srf_generation.realisation import RealisationFault

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
UTM_CRS = pyproj.CRS(proj="utm", zone=11, ellps="WGS84")
WGS2UTM = pyproj.Transformer.from_crs(WGS_CODE, UTM_CRS)
TRANSFORMER_MAP = {"nztm": WGS2NZTM, "utm": WGS2UTM}
# Convert x, y to lat, lon
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)
SUBDIVISION_RESOLUTION_KM = 0.1
FAULTSEG2GSFDIPDIR = "fault_seg2gsf_dipdir"
SRF2STOCH = "srf2stoch"


def normalise_name(name: str) -> str:
    """Normalise a fault name for yaml file output.

    Parameters
    ----------
    name : str
        The name to normalise.

    Returns
    -------
    str
        The normalised version of the name.
    """
    fault_no_illegal_characters = re.sub(r" |,|:", "_", name)
    return fault_no_illegal_characters.lower()


def generate_fault_gsf(gsf_output_directory: Path, fault: RealisationFault):
    gsf_output_filepath = gsf_output_directory / f"{normalise_name(fault.name)}.gsf"
    gsf_df = pd.DataFrame(
        [
            {
                "strike": plane.strike,
                "dip": plane.dip,
                "length": plane.length,
                "width": plane.width,
                "rake": plane.rake,
                "meshgrid": gsf.coordinate_meshgrid(
                    plane.corners[0],
                    plane.corners[1],
                    plane.corners[-1],
                    SUBDIVISION_RESOLUTION_KM * 1000,
                ),
            }
            for plane in fault.planes
        ]
    )
    gsf.write_fault_to_gsf_file(
        gsf_output_filepath, gsf_df, SUBDIVISION_RESOLUTION_KM * 1000
    )

    return gsf_output_filepath


def srf_file_for_fault(output_directory: Path, fault: RealisationFault) -> Path:
    return output_directory / (fault.name + ".srf")


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


def generate_type4_fault_srf(
    realisation: realisation.Realisation,
    fault: RealisationFault,
    output_directory: Path,
):
    gsf_output_directory = output_directory / "gsf"
    gsf_file_path = generate_fault_gsf(gsf_output_directory, fault)

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)

    resolution = 100
    nx = sum(
        gsf.gridpoint_count_in_length(plane.length_m, resolution)
        for plane in fault.planes
    )
    ny = gsf.gridpoint_count_in_length(fault.planes[0].width_m, resolution)
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
    print(" ".join(genslip_cmd))
    srf_file_path = srf_file_for_fault(output_directory, fault)
    with open(srf_file_path, "w", encoding="utf-8") as srf_file_handle:
        subprocess.run(
            genslip_cmd, stdout=srf_file_handle, stderr=subprocess.PIPE, check=True
        )


def closest_gridpoint(grid: np.ndarray, point: np.ndarray) -> int:
    return int(np.argmin(np.sum(np.square(grid - point), axis=1)))


def stitch_srf_files(
    realisation_obj: realisation.Realisation, output_directory: Path
) -> Path:
    srf_output_filepath = output_directory / f"{realisation_obj.name}.srf"
    with open(srf_output_filepath, "w", encoding="utf-8") as srf_file_output:
        fault_points = {}
        header = []

        srf.write_version(srf_file_output)

        for fault in realisation.topologically_sorted_faults(
            realisation_obj.faults.values()
        ):
            with open(
                srf_file_for_fault(output_directory, fault), "r", encoding="utf-8"
            ) as fault_srf_file:
                srf.read_version(fault_srf_file)
                fault_header = srf.read_srf_headers(fault_srf_file)
                if fault.parent:
                    for plane in fault_header:
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
                parent_coords = wgsdepth_to_nztm(
                    np.array(fault.parent_jump_coords).reshape((1, -1))
                ).ravel()
                parent_fault_points = fault_points[parent.name]
                grid_points = wgsdepth_to_nztm(
                    np.array(
                        [
                            [point.lat, point.lon, point.dep * 1000]
                            for point in parent_fault_points
                        ]
                    )
                )
                jump_index = closest_gridpoint(
                    grid_points,
                    parent_coords,
                )
                t_delay = parent_fault_points[jump_index].tinit
            cur_fault_points = fault_points[fault.name]
            for point in cur_fault_points:
                point.tinit += t_delay
                srf.write_srf_point(srf_file_output, point)

        return srf_output_filepath


def generate_type4_fault_srfs_parallel(
    realisation: realisation.Realisation,
    output_directory: Path,
):
    # need to do this before multiprocessing because of race conditions
    gsf_directory = output_directory / "gsf"
    if not gsf_directory.exists():
        gsf_directory.mkdir()
    with multiprocessing.Pool() as worker_pool:
        worker_pool.starmap(
            generate_type4_fault_srf,
            [
                (
                    realisation,
                    fault,
                    output_directory,
                )
                for fault in realisation.faults.values()
            ],
        )


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
):
    """Generate a type-5 SRF file from a given realisation specification."""
    realisation_parameters = realisation.read_realisation(realisation_filepath)
    generate_type4_fault_srfs_parallel(realisation_parameters, output_directory)
    srf_output_filepath = stitch_srf_files(realisation_parameters, output_directory)
    create_stoch(
        srf_output_filepath, output_directory / f"{realisation_parameters.name}.stoch"
    )


if __name__ == "__main__":
    typer.run(main)

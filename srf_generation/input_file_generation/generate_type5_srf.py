import argparse
import collections
import multiprocessing
import subprocess
import tempfile
from pathlib import Path
from typing import Annotated, Generator

import numpy as np
import pandas as pd
import pyproj
import typer
from qcore import binary_version
from srf_generation import realisation
from srf_generation.realisation import RealisationFault

import srf

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


def type4_fault_gsf(fault: RealisationFault):
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as gsf_output_file:
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
            gsf_output_file.name, gsf_df, SUBDIVISION_RESOLUTION_KM * 1000
        )
        print(gsf_output_file.name)

    return gsf_output_file.name


def srf_file_for_fault(output_directory: Path, fault: RealisationFault) -> Path:
    return output_directory / (fault.name + ".srf")


def create_stoch(
    stoch_file: str,
    srf_file: str,
):
    """Create a stoch file from a SRF file.

    Parameters
    ----------
    stoch_file : str
        The filepath to output the stoch file to.
    srf_file : str
        The filepath of the SRF file.
    single_segment : bool
        True if the stoch file is a single segment.
    """

    dx = 2.0
    dy = 2.0
    srf2stoch = binary_version.get_unversioned_bin(SRF2STOCH)
    command = [
        srf2stoch,
        f"dx={dx}",
        f"dy={dy}",
        f"infile={srf_file}",
        f"outfile={stoch_file}",
    ]
    subprocess.run(command, stderr=subprocess.PIPE, check=True)


def generate_type4_fault_srf(
    realisation: realisation.Realisation,
    fault: RealisationFault,
    output_directory: Path,
):
    gsf_file_path = type4_fault_gsf(fault)

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)

    resolution = 100
    nx = sum(
        gsf.gridpoints_along_direction(plane.length_m, resolution)
        for plane in fault.planes
    )
    ny = gsf.gridpoints_along_direction(fault.planes[0].width_m, resolution)
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


def topologically_sorted_faults(
    faults: list[RealisationFault],
) -> Generator[RealisationFault, None, None]:
    fault_children_map = collections.defaultdict(list)
    for fault in faults:
        if fault.parent:
            fault_children_map[fault.parent.name].append(fault)

    def bfs_traverse_fault_map(
        fault: RealisationFault,
    ) -> Generator[RealisationFault, None, None]:
        yield fault
        for child in fault_children_map[fault.name]:
            yield from bfs_traverse_fault_map(child)

    initial_fault = next(fault for fault in faults if not fault.parent)
    yield from bfs_traverse_fault_map(initial_fault)


def stitch_srf_files(realisation: realisation.Realisation, output_directory: Path):
    srf_output_filepath = output_directory / f"{realisation.name}.srf"
    with open(srf_output_filepath, "w") as srf_file_output:
        fault_srfs = {}
        fault_points = {}
        header = []

        srf.write_version(srf_file_output)

        for fault in topologically_sorted_faults(realisation.faults.values()):
            fault_srf_file = open(srf_file_for_fault(output_directory, fault), "r")
            fault_srfs[fault.name] = fault_srf_file
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
        for fault in topologically_sorted_faults(realisation.faults.values()):
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
                print(
                    f"Parent point lat/lon",
                    nztm_to_wgsdepth(parent_coords.reshape((1, -1))),
                )
                print(
                    f"Closest parent gridpoint: {nztm_to_wgsdepth(grid_points[jump_index].reshape((1, -1)))}"
                )
                t_delay = parent_fault_points[jump_index].tinit
                print(f"Fault {fault.name} will have a delay of {t_delay}")
            cur_fault_points = fault_points[fault.name]
            for point in cur_fault_points:
                point.tinit += t_delay
                srf.write_srf_point(srf_file_output, point)

        for fault_srf_file in fault_srfs.values():
            fault_srf_file.close()


def generate_type4_fault_srfs_parallel(
    realisation: realisation.Realisation,
    output_directory: Path,
):
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


def wgsdepth_to_nztm(wgsdepthcoordinates: np.ndarray) -> np.ndarray:
    nztm_coords = np.array(
        WGS2NZTM.transform(wgsdepthcoordinates[:, 0], wgsdepthcoordinates[:, 1]),
    ).T
    return np.append(nztm_coords, wgsdepthcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def nztm_to_wgsdepth(nztmcoordinates: np.ndarray) -> np.ndarray:
    wgs_coords = np.array(
        NZTM2WGS.transform(nztmcoordinates[:, 0], nztmcoordinates[:, 1]),
    ).T
    return np.append(wgs_coords, nztmcoordinates[:, 2].reshape((-1, 1)), axis=-1)


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
    stitch_srf_files(realisation, output_directory)


if __name__ == "__main__":
    typer.run(main)

import dataclasses
import multiprocessing
import subprocess
import tempfile
from pathlib import Path
import fault
import pyproj
import pdb

import numpy as np
import qcore.geo
import scipy as sp
import yaml
from qcore import binary_version
from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)

import rupture_propogation
from rupture_propogation import RuptureCausalityTree

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


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    magnitude: float
    # moment: float
    dt: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    initial_fault: str
    shypo: float
    dhypo: float
    velocity_model: str
    faults: dict[str, fault.Fault]

    def rupture_area(self) -> float:
        return sum(fault.area() for fault in self.faults.values())

    def fault_magnitude_by_area(self, fault: str) -> float:
        return self.magnitude * (self.faults[fault].area() / self.rupture_area())


def read_realisation(realisation_filepath: Path) -> Realisation:
    with open(realisation_filepath, "r", encoding="utf-8") as realisation_file:
        raw_yaml_data = yaml.safe_load(realisation_file)
        faults = {
            name: fault.Fault(
                name=name,
                tect_type=fault_obj["tect_type"],
                segments=[
                    fault.FaultSegment(**params) for params in fault_obj["segments"]
                ],
            )
            for name, fault_obj in raw_yaml_data.pop("faults").items()
        }
        return Realisation(
            **raw_yaml_data,
            faults=faults,
            velocity_model=DEFAULT_1D_VELOCITY_MODEL_PATH,
        )


def type4_fault_gsf(realisation: Realisation, fault: fault.Fault):
    fault_seg_bin = binary_version.get_unversioned_bin(FAULTSEG2GSFDIPDIR)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as input_file, tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as gsf_output_file:
        input_file.write(f"{len(fault.segments)}\n")
        for segment in fault.segments:
            input_file.write(
                f"{segment.clon:6f} {segment.clat:6f} {segment.dtop:6f} {segment.strike:.4f} {segment.dip:.4f} {segment.rake:.4f} {segment.length:.4f} {segment.width:.4f} {segment.length_subdivisions():d} {segment.width_subdivisions():d}\n"
            )
        input_file.flush()
        fault_seg_command = [
            fault_seg_bin,
            "read_slip_vals=0",
            f"infile={input_file.name}",
            f"outfile={gsf_output_file.name}",
            f"dipdir={segment.dip_dir}",
        ]
        print(" ".join(fault_seg_command))
        subprocess.run(fault_seg_command, check=True)

    return gsf_output_file.name


def generate_type4_fault_srf(
    realisation: Realisation,
    fault_name: str,
    output_directory: Path,
    hypocentre: np.ndarray,
):
    fault = realisation.faults[fault_name]

    gsf_file_path = type4_fault_gsf(realisation, fault)

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)
    magnitude = realisation.fault_magnitude_by_area(fault_name)

    lengths = fault.lengths()
    nx = int(np.sum(np.round(lengths / SUBDIVISION_RESOLUTION_KM)))
    ny = int(np.round(fault.widths()[0] / SUBDIVISION_RESOLUTION_KM))

    genslip_cmd = [
        genslip_bin,
        "read_erf=0",
        "write_srf=1",
        "read_gsf=1",
        "write_gsf=0",
        f"infile={gsf_file_path}",
        f"mag={magnitude}",
        f"nstk={nx}",
        f"ndip={ny}",
        "ns=1",
        "nh=1",
        f"seed={realisation.genslip_seed}",
        f"velfile={realisation.velocity_model}",
        f"shypo={hypocentre[0]}",
        f"dhypo={hypocentre[1]}",
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

    srf_file_path = output_directory / (fault_name + ".srf")
    with open(srf_file_path, "w", encoding="utf-8") as srf_file_handle:
        subprocess.run(
            genslip_cmd, stdout=srf_file_handle, stderr=subprocess.PIPE, check=True
        )


def generate_type4_fault_srfs_parallel(
    realisation: Realisation,
    hypocentres: dict[str, (float, float)],
    output_directory: Path,
):
    with multiprocessing.Pool() as worker_pool:
        worker_pool.starmap(
            generate_type4_fault_srf,
            [
                (realisation, fault_name, output_directory, hypocentres[fault_name])
                for fault_name in realisation.faults
            ],
        )


def wgsdepth_to_nztm(wgsdepthcoordinates: np.ndarray) -> np.ndarray:
    nztm_coords = np.array(
        WGS2NZTM.transform(wgsdepthcoordinates[:, 0], wgsdepthcoordinates[:, 1]),
    ).T
    return np.append(nztm_coords, wgsdepthcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def closest_points_between_faults(fault_u, fault_v):
    scaling_factor = np.array([1e3, 1e3, 1])
    fault_u_corners = (
        np.array([wgsdepth_to_nztm(corner) for corner in fault_u.corners()])
        / scaling_factor
    )
    fault_v_corners = (
        np.array([wgsdepth_to_nztm(corner) for corner in fault_v.corners()])
        / scaling_factor
    )
    point_u, point_v = qcore.geo.closest_points_between_plane_sequences(
        fault_u_corners, fault_v_corners
    )
    return point_u * scaling_factor, point_v * scaling_factor


def build_rupture_causality_tree(realisition: Realisation) -> RuptureCausalityTree:
    jump_point_map = {
        fault_u_name: {
            fault_v_name: closest_points_between_faults(fault_u, fault_v)
            for fault_v_name, fault_v in realisation.faults.items()
            if fault_v_name != fault_u_name
        }
        for fault_u_name, fault_u in realisation.faults.items()
    }
    distance_graph = {
        fault_u_name: {
            fault_v_name: sp.spatial.distance.cdist(
                u_point.reshape((1, -1)), v_point.reshape((1, -1))
            )[0, 0]
            / 1000
            for fault_v_name, (u_point, v_point) in jump_point_map[fault_u_name].items()
        }
        for fault_u_name in jump_point_map
    }
    print(distance_graph)

    # pruned = rupture_propogation.prune_distance_graph(distance_graph, 15)
    probability_graph = rupture_propogation.probability_graph(distance_graph)

    return rupture_propogation.probabilistic_shortest_path(
        probability_graph, realisation.initial_fault
    )

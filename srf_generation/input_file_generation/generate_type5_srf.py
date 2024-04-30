import dataclasses
import multiprocessing
import subprocess
import tempfile
from pathlib import Path
import fault
import collections
import pyproj
import itertools


import numpy as np
import qcore.geo
import scipy as sp
import yaml
import srf
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
class FaultJump:
    parent: fault.Fault
    jump_location_lat: float
    jump_location_lon: float


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
    velocity_model: str
    faults: dict[str, fault.Fault]
    causality_map: dict[str, FaultJump]

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
                shyp=fault_obj["shyp"],
                dhyp=fault_obj["dhyp"],
            )
            for name, fault_obj in raw_yaml_data.pop("faults").items()
        }

        causality_map = {}
        for fault_name, jump_obj in raw_yaml_data["causality_map"].items():
            if jump_obj["parent"] is None:
                causality_map[fault_name] = None
            else:
                causality_map[fault_name] = FaultJump(**jump_obj)

        return Realisation(
            name=raw_yaml_data["name"],
            type=raw_yaml_data["type"],
            dt=raw_yaml_data["dt"],
            genslip_seed=raw_yaml_data["genslip_seed"],
            srfgen_seed=raw_yaml_data["srfgen_seed"],
            initial_fault=raw_yaml_data["initial_fault"],
            genslip_version=raw_yaml_data["genslip_version"],
            magnitude=raw_yaml_data["magnitude"],
            faults=faults,
            velocity_model=DEFAULT_1D_VELOCITY_MODEL_PATH,
            causality_map=causality_map,
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
        subprocess.run(fault_seg_command, check=True)

    return gsf_output_file.name


def srf_file_by_name(output_directory: Path, fault_name: str) -> Path:
    return output_directory / (fault_name + ".srf")


def generate_type4_fault_srf(
    realisation: Realisation,
    fault_name: str,
    output_directory: Path,
):
    fault = realisation.faults[fault_name]

    gsf_file_path = type4_fault_gsf(realisation, fault)

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)
    magnitude = realisation.fault_magnitude_by_area(fault_name)
    print(f"{fault.name}: {magnitude}")

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
        f"mag={realisation.magnitude}",
        f"nstk={nx}",
        f"ndip={ny}",
        "ns=1",
        "nh=1",
        f"seed={realisation.genslip_seed}",
        f"velfile={realisation.velocity_model}",
        f"shypo={fault.shyp}",
        f"dhypo={fault.dhyp}",
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
    srf_file_path = srf_file_by_name(output_directory, fault.name)
    with open(srf_file_path, "w", encoding="utf-8") as srf_file_handle:
        subprocess.run(
            genslip_cmd, stdout=srf_file_handle, stderr=subprocess.PIPE, check=True
        )


def stitch_srf_files(realisation: Realisation, output_directory: Path):
    srf_output_filepath = srf_file_by_name(output_directory, realisation.name)
    with open(srf_output_filepath, "w") as srf_file_output:
        fault_srfs = {}
        fault_points = {}
        header = []

        srf.write_version(srf_file_output)

        for fault_name in realisation.faults:
            fault_srf_file = open(srf_file_by_name(output_directory, fault_name), "r")
            fault_srfs[fault_name] = fault_srf_file
            srf.read_version(fault_srf_file)
            fault_header = srf.read_srf_headers(fault_srf_file)
            if fault_name != realisation.initial_fault:
                for segment in fault_header:
                    segment.shyp = -999
                    segment.dhyp = -999
            header.extend(fault_header)
            point_count = srf.read_points_count(fault_srf_file)
            fault_points[fault_name] = srf.read_srf_n_points(
                point_count, fault_srf_file
            )

        srf.write_srf_header(srf_file_output, header)
        srf.write_point_count(
            srf_file_output, sum(len(points) for points in fault_points.values())
        )

        for fault_name, fault_points in fault_points.items():
            t_delay = 0
            if fault_name != realisation.initial_fault:
                # find closest grid point to the jump location
                # compute the time delay as equal to the tinit of this point (for now)
                pass
            for point in fault_points:
                point.tinit += t_delay
                srf.write_srf_point(srf_file_output, point)

        for fault_srf_file in fault_srfs.values():
            fault_srf_file.close()


def generate_type4_fault_srfs_parallel(
    realisation: Realisation,
    output_directory: Path,
):
    with multiprocessing.Pool() as worker_pool:
        worker_pool.starmap(
            generate_type4_fault_srf,
            [
                (
                    realisation,
                    fault_name,
                    output_directory,
                )
                for fault_name in realisation.faults
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


def closest_points_between_faults(fault_u, fault_v):
    fault_u_corners = np.array(
        [wgsdepth_to_nztm(corner) for corner in fault_u.corners()]
    )

    fault_v_corners = np.array(
        [wgsdepth_to_nztm(corner) for corner in fault_v.corners()]
    )

    point_u, point_v = qcore.geo.closest_points_between_plane_sequences(
        fault_u_corners, fault_v_corners
    )
    return point_u, point_v


def test_fault_segment_vialibility(
    fault1: fault.FaultSegment, fault2: fault.FaultSegment, cutoff: float
) -> bool:
    fault1_centroid = wgsdepth_to_nztm(fault1.centroid().reshape((1, -1))).ravel()
    fault2_centroid = wgsdepth_to_nztm(fault2.centroid().reshape((1, -1))).ravel()
    fault1_radius = max(fault1.width_m, fault1.length_m) + cutoff
    fault2_radius = max(fault2.width_m, fault2.length_m) + cutoff
    return qcore.geo.spheres_intersect(
        fault1_centroid, fault1_radius, fault2_centroid, fault2_radius
    )


def test_fault_viability(
    fault1: fault.Fault, fault2: fault.Fault, cutoff: float
) -> bool:
    return any(
        test_fault_segment_vialibility(seg1, seg2, cutoff)
        for (seg1, seg2) in itertools.product(fault1.segments, fault2.segments)
    )


def build_rupture_causality_tree(realisation: Realisation) -> RuptureCausalityTree:
    jump_point_map = collections.defaultdict(dict)
    cutoff = 15000
    for fault_u_name, fault_u in realisation.faults.items():
        for fault_v_name, fault_v in realisation.faults.items():
            if fault_u_name == fault_v_name:
                continue
            if test_fault_viability(fault_u, fault_v, cutoff):
                jump_point_map[fault_u_name][fault_v_name] = (
                    closest_points_between_faults(fault_u, fault_v)
                )

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

    pruned = rupture_propogation.prune_distance_graph(distance_graph, 15000)
    probability_graph = rupture_propogation.probability_graph(pruned)

    return rupture_propogation.probabilistic_minimum_spanning_tree(
        probability_graph, realisation.initial_fault
    )


def compute_jump_point_hypocentre_fault_coordinates(
    from_fault: fault.Fault, to_fault: fault.Fault
) -> (float, float):
    _, to_fault_point_nztm = closest_points_between_faults(from_fault, to_fault)
    return to_fault.hypocentre_wgs_to_fault_coordinates(
        nztm_to_wgsdepth(to_fault_point_nztm.reshape((1, -1))).ravel()
    )


def compute_hypocentres(
    realisation: Realisation, rupture_causality_tree: RuptureCausalityTree
):
    hypocentres = {}
    for fault_name, to_fault in realisation.faults.items():
        if rupture_causality_tree[fault_name] is None:
            hypocentres[fault_name] = (realisation.shypo, realisation.dhypo)
        else:
            from_fault = realisation.faults[rupture_causality_tree[fault_name]]
            hypocentres[fault_name] = compute_jump_point_hypocentre_fault_coordinates(
                from_fault, to_fault
            )
    return hypocentres


if __name__ == "__main__":
    realisation_filepath = "/home/jake/Downloads/good_rupture.yaml"
    realisation = read_realisation(realisation_filepath)
    rupture_causality_tree = build_rupture_causality_tree(realisation)
    hypocentres = compute_hypocentres(realisation, rupture_causality_tree)
    print(hypocentres)

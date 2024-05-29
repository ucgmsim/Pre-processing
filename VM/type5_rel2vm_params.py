#!/usr/bin/env python3


# Initial screening:
# Need at least min_points_above_thresh points in the interval
# initial_screening_min_freq_Hz to initial_screening_max_freq_Hz with
# SNR > initial_screening_snr_thresh_ver for the ver component and
# SNR > initial_screening_snr_thresh_horiz for the 000 and 090 components
#
#
#
def find_fmax(
    filename: Path, metadata: pd.DataFrame
) -> tuple[dict[str, Any], dict[str, Any]]:
    """some short description.

    Parameters
    ----------
    filename : Path
        Path to the ___snr_fas.csv file.
    metadata : pd.DataFrame
        Contains the SNR meta data.

    Returns
    -------
    fmax_record : dict[str, Any] or None
        Dictionary with keys record_id, fmax_000, fmax_090, fmax_ver. This value
        will be None when a record is skipped or not recorded.
    skipped_record : dict[str, Any] or None
        Records a skipped record in a dictionary with keys record_id and reason,
        if a skipped record exists.
    """
    pass


import dataclasses
from pathlib import Path
from typing import Annotated, Tuple

import numpy as np
import scipy as sp
import typer
import yaml
from qcore import coordinates, geo
from srf_generation import realisation

from models.AfshariStewart_2016_Ds import Afshari_Stewart_2016_Ds
from models.classdef import Fault, FaultStyle, Site, TectType, estimate_z1p0
from VM.models.Bradley_2010_Sa import Bradley_2010_Sa


@dataclasses.dataclass
class BoundingBox:
    corners: np.ndarray

    @property
    def origin(self):
        return np.mean(self.corners, axis=0)

    @property
    def extent_x(self):
        return np.linalg.norm(np.linalg.norm(self.corners[2] - self.corners[1]) / 1000)

    @property
    def extent_y(self):
        return np.linalg.norm(np.linalg.norm(self.corners[1] - self.corners[0]) / 1000)

    @property
    def bearing(self):

        north_direction = np.array([1, 0, 0])
        up_direction = np.array([0, 0, 1])
        horizontal_direction = np.append(self.corners[1] - self.corners[0], 0)
        return geo.oriented_bearing_wrt_normal(
            north_direction, horizontal_direction, up_direction
        )


def axis_aligned_bounding_box(points: np.ndarray) -> np.ndarray:
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    return np.array([[min_x, min_y], [max_x, min_y], [max_x, max_y], [min_x, max_y]])


def rotation_matrix(angle: float) -> np.ndarray:
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


def minimum_area_bounding_box(points: np.ndarray) -> BoundingBox:
    """Find the minimum area bounding rectangle surrounding the points."""
    convex_hull = sp.spatial.ConvexHull(points).points
    segments = np.array(
        [
            convex_hull[(i + 1) % len(convex_hull)] - convex_hull[i]
            for i in range(len(convex_hull))
        ]
    )
    rotation_angles = -np.arctan2(segments[:, 1], segments[:, 0])

    bounding_boxes = [
        axis_aligned_bounding_box(convex_hull @ rotation_matrix(angle).T)
        for angle in rotation_angles
    ]

    rotation_angle, minimum_bounding_box = min(
        zip(rotation_angles, bounding_boxes),
        key=lambda box: np.linalg.norm(box[1][1] - box[1][0])
        * np.linalg.norm(box[1][2] - box[1][1]),
    )
    return BoundingBox(minimum_bounding_box @ rotation_matrix(-rotation_angle).T)


def realisation_bounding_box(type5_realisation: realisation.Realisation) -> BoundingBox:
    fault_corners = np.vstack(
        [fault.corners_nztm() for fault in type5_realisation.faults.values()]
    )[:, :2]
    return minimum_area_bounding_box(fault_corners)


def mag2pgv(mag: float):
    """
    Return PGV for a given magnitude based on default scaling relationship
    Parameters
    ----------
    mag : magnitude

    Returns
    -------
    Interpolated pgv

    """
    return np.interp(
        mag,
        [3.5, 4.1, 4.7, 5.2, 5.5, 5.8, 6.2, 6.5, 6.8, 7.0, 7.4, 7.7, 8.0],
        [0.015, 0.0375, 0.075, 0.15, 0.25, 0.4, 0.7, 1.0, 1.35, 1.65, 2.1, 2.5, 3.0],
    )


def find_rrup(magnitude: float) -> Tuple[float, float]:
    """
    Find rrup at which pgv estimated from magnitude is close to target.

    Parameters
    ----------
    pgv_target : pgv that we wish to come close to

    Returns
    -------

    """

    pgv_target = mag2pgv(magnitude)

    def pgv_delta_from_rrup(rrup: float):
        siteprop = Site()
        siteprop.vs30 = 500
        siteprop.Rtvz = 0
        siteprop.vs30measured = False
        siteprop.z1p0 = estimate_z1p0(siteprop.vs30)
        faultprop = Fault()
        faultprop.dip = 0
        faultprop.rake = 0
        faultprop.Mw = magnitude
        faultprop.ztor = (
            0.0  # DON'T CHANGE THIS - this assumes we have a surface point source
        )
        faultprop.tect_type = TectType.ACTIVE_SHALLOW
        faultprop.faultstyle = FaultStyle.UNKNOWN
        siteprop.Rrup = rrup
        siteprop.Rx = rrup
        siteprop.Rjb = rrup
        return np.abs(Bradley_2010_Sa(siteprop, faultprop, "PGV")[0] - pgv_target)

    rrup_optimise_result = sp.optimize.minimize_scalar(
        pgv_delta_from_rrup, bounds=(0, 1e4)
    )
    rrup = rrup_optimise_result.x
    pgv_delta = rrup_optimise_result.fun
    return rrup, pgv_target + pgv_delta


def grow_bounding_box(bounding_box: BoundingBox, radius: float):
    """Grow a realisation bounding box so that bounding box contains a circle of
    given radius centred on the origin."""
    centroid = bounding_box.origin
    # basically just guessing at this point
    corners_centred = bounding_box.corners - centroid
    centre_x = (corners_centred[2] - corners_centred[1]) / 2
    centre_y = (corners_centred[1] - corners_centred[0]) / 2
    scale_x = radius / (np.linalg.norm(centre_x) / 1000)
    scale_y = radius / (np.linalg.norm(centre_y) / 1000)

    scale = min(scale_x, scale_y)
    basis = np.array([centre_x, centre_y])
    return BoundingBox(
        scale
        * np.array(
            [
                [-1, -1],
                [-1, 1],
                [1, 1],
                [1, -1],
            ]
        )
        @ basis
        + centroid
    )


def guess_simulation_duration(
    bounding_box: BoundingBox,
    type5_realisation: realisation.Realisation,
    ds_multiplier: float,
) -> float:
    """Guess how long a simulation should go based on bounding box of the domain and the location of the faults."""
    fault_corners = np.vstack(
        [fault.corners_nztm() for fault in type5_realisation.faults.values()]
    )
    fault_centroid = np.mean(fault_corners, axis=0)
    box_corners = np.append(
        bounding_box.corners,
        np.zeros((4, 1)),
        axis=1,
    )
    maximum_distance = np.max(np.linalg.norm(box_corners - fault_centroid, axis=1))
    s_wave_m_per_s = 3500
    s_wave_arrival_time = maximum_distance / s_wave_m_per_s

    # compute the pairwise distance between the domain corners and the fault corners
    pairwise_distance = np.linalg.norm(
        box_corners[np.newaxis, :, :] - fault_corners[:, np.newaxis, :], axis=2
    )
    largest_corner_distance = np.max(np.min(pairwise_distance, axis=1)) / 1000

    # taken from rel2vm_params.py
    siteprop = Site()
    siteprop.vs30 = 500
    siteprop.Rtvz = 0
    siteprop.vs30measured = False
    siteprop.z1p0 = estimate_z1p0(siteprop.vs30)
    siteprop.Rrup = largest_corner_distance
    faultprop = Fault()
    faultprop.ztor = (
        0.0  # DON'T CHANGE THIS - this assumes we have a surface point source
    )
    faultprop.tect_type = TectType.ACTIVE_SHALLOW
    faultprop.faultstyle = FaultStyle.UNKNOWN
    faultprop.Mw = type5_realisation.magnitude

    ds = Afshari_Stewart_2016_Ds(siteprop, faultprop, "Ds595")[0]

    return s_wave_arrival_time + ds_multiplier * ds


def get_max_depth(mag: float, depth: float):
    """
    Maximum depth (ie. z extent) based on magnitude, hypocentre depth.

    Parameters
    ----------
    mag : magnitude
    depth : hypocentre depth

    Returns
    -------
    maximum depth
    """
    return round(
        10
        + depth
        + (10 * np.power((0.5 * np.power(10, (0.55 * mag - 1.2)) / depth), 0.3)),
        0,
    )


def get_rupture_corners(db, rupture_id: int) -> np.ndarray:
    return np.vstack(
        [fault.corners_nztm() for fault in db.get_rupture_faults(rupture_id)]
    )[:, :2]


def rupture_area(db, rupture_id: int) -> float:
    bounding_box = minimum_area_bounding_box(get_rupture_corners(db, rupture_id))

    return max(bounding_box.extent_x, bounding_box.extent_y)


def main(
    realisation_filepath: Annotated[
        Path,
        typer.Argument(
            help="The path to the realisation to generate VM parameters for."
        ),
    ],
    output_path: Annotated[
        Path, typer.Argument(help="The path to the output VM params.")
    ],
    resolution: Annotated[
        float, typer.Option(help="The resolution of the simulation in kilometres.")
    ] = 0.1,
):
    "Generate velocity model parameters for a Type-5 realisation."
    type5_realisation = realisation.read_realisation(realisation_filepath)
    bounding_box = realisation_bounding_box(type5_realisation)
    min_depth = 0
    initial_fault = type5_realisation.initial_fault()
    magnitude = type5_realisation.magnitude
    max_depth = get_max_depth(
        magnitude,
        initial_fault.fault_coordinates_to_wgsdepth_coordinates(
            np.array([initial_fault.shyp, initial_fault.dhyp])
        )[2]
        / 1000,
    )
    nx = int(np.ceil(bounding_box.extent_x / resolution))
    ny = int(np.ceil(bounding_box.extent_y / resolution))
    nz = int(np.ceil((max_depth - min_depth) / (resolution)))
    sim_duration = guess_simulation_duration(bounding_box, type5_realisation, 1.2)
    (rrup, _) = find_rrup(type5_realisation.magnitude)
    bounding_box = grow_bounding_box(bounding_box, rrup)
    box_origin_coords = coordinates.nztm_to_wgs_depth(np.append(bounding_box.origin, 0))
    vm_params = {
        "mag": magnitude,
        "MODEL_LAT": float(box_origin_coords[0]),
        "MODEL_LON": float(box_origin_coords[1]),
        "MODEL_ROT": float(bounding_box.bearing),
        "hh": resolution,
        "extent_x": float(bounding_box.extent_x),
        "extent_y": float(bounding_box.extent_y),
        "extent_zmax": float(max_depth),
        "extent_zmin": float(min_depth),
        "sim_duration": float(sim_duration),
        "nx": nx,
        "ny": ny,
        "nz": nz,
    }
    with open(output_path, "w") as output_file_handle:
        yaml.dump(vm_params, output_file_handle)


if __name__ == "__main__":
    typer.run(main)

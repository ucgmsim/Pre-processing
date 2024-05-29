#!/usr/bin/env python3

import dataclasses
from pathlib import Path
from typing import Annotated

import numpy as np
import scipy as sp
import typer
import yaml
from qcore import coordinates, geo
from srf_generation import realisation

from VM.models.AfshariStewart_2016_Ds import Afshari_Stewart_2016_Ds
from VM.models.classdef import Fault, FaultStyle, Site, TectType, estimate_z1p0


def axis_aligned_bounding_box(points: np.ndarray) -> np.ndarray:
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    return np.array([[min_x, min_y], [max_x, min_y], [max_x, max_y], [min_x, max_y]])


def rotation_matrix(angle: float) -> np.ndarray:
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


def minimum_area_bounding_box(points: np.ndarray) -> np.ndarray:
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
    return minimum_bounding_box @ rotation_matrix(-rotation_angle).T


@dataclasses.dataclass
class BoundingBox:
    origin: np.ndarray
    bearing: float
    extent_x: float
    extent_y: float
    corners: np.ndarray


def bounding_box_bearing(rectangle: np.ndarray) -> float:
    """Find the bearing of a bounding box from north"""
    north_direction = np.array([1, 0, 0])
    up_direction = np.array([0, 0, 1])
    horizontal_direction = rectangle[1] - rectangle[0]
    return geo.oriented_bearing_wrt_normal(
        north_direction, horizontal_direction, up_direction
    )


def realisation_bounding_box(type5_realisation: realisation.Realisation) -> BoundingBox:
    fault_corners = np.vstack(
        [fault.corners_nztm() for fault in type5_realisation.faults.values()]
    )
    bounding_box_corners = minimum_area_bounding_box(fault_corners)
    bearing = bounding_box_bearing(bounding_box_corners)
    origin = np.mean(bounding_box_corners)
    extent_x = np.linalg.norm(bounding_box_corners[1] - bounding_box_corners[0]) / 1000
    extent_y = np.linalg.norm(bounding_box_corners[2] - bounding_box_corners[1]) / 1000
    return BoundingBox(origin, bearing, extent_x, extent_y, bounding_box_corners)


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
        bounding_box.corners, np.zeros_like(bounding_box.corners[:, 0])
    )
    maximum_distance = np.max(np.linalg.norm(box_corners - fault_centroid, axis=1))
    s_wave_km_per_s = 3.5
    s_wave_arrival_time = maximum_distance / s_wave_km_per_s

    # compute the pairwise distance between the domain corners and the fault corners
    pairwise_distance = np.linalg.norm(
        box_corners[np.newaxis, :, :] - fault_corners[:, np.newaxis, :], axis=2
    )
    largest_corner_distance = np.max(np.min(pairwise_distance, axis=1))

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
        initial_fault.fault_coordinates_to_global_coordinates(
            np.array([initial_fault.shyp, initial_fault.dhyp])
        ),
    )
    nx = int(np.ceil(bounding_box.extent_x / resolution))
    ny = int(np.ceil(bounding_box.extent_y / resolution))
    nz = int(np.ceil((max_depth - min_depth) / resolution))
    sim_duration = guess_simulation_duration(type5_realisation)
    box_origin_coords = coordinates.nztm_to_wgs_depth(np.append(bounding_box.origin, 0))
    vm_params = {
        "mag": magnitude,
        "MODEL_LAT": box_origin_coords[0],
        "MODEL_LON": box_origin_coords[1],
        "MODEL_ROT": bounding_box.bearing,
        "hh": resolution,
        "extent_x": bounding_box.extent_x,
        "extent_y": bounding_box.extent_y,
        "extent_zmax": max_depth,
        "extent_zmin": min_depth,
        "sim_duration": sim_duration,
        "nx": nx,
        "ny": ny,
        "nz": nz,
    }
    with open(output_path, "w") as output_file_handle:
        yaml.dump(vm_params, output_file_handle)


if __name__ == "__main__":
    typer.run(main)

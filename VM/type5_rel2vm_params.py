#!/usr/bin/env python3

import dataclasses
from pathlib import Path
from typing import Annotated, Tuple

import numpy as np
import scipy as sp
import shapely
import typer
import yaml
from models.AfshariStewart_2016_Ds import Afshari_Stewart_2016_Ds
from models.classdef import Fault, FaultStyle, Site, TectType, estimate_z1p0
from qcore import coordinates, geo
from shapely import Polygon

from srf_generation import realisation
from VM import bounding_box
from VM.bounding_box import BoundingBox
from VM.models.Bradley_2010_Sa import Bradley_2010_Sa

script_dir = Path(__file__).resolve().parent
NZ_LAND_OUTLINE = script_dir / "../SrfGen/NHM/res/rough_land.txt"


def get_land_polygon():
    polygon_coordinates = np.loadtxt(NZ_LAND_OUTLINE)
    polygon_coordinates = polygon_coordinates[:, ::-1]
    polygon_coordinates = np.append(
        polygon_coordinates,
        np.zeros_like(polygon_coordinates[:, 0]).reshape((-1, 1)),
        axis=1,
    )
    polygon_coordinates = coordinates.wgs_depth_to_nztm(polygon_coordinates)[:, :2]
    return shapely.Polygon(polygon_coordinates)


def pgv_estimate_from_magnitude(magnitude: float):
    """
    Return PGV for a given magnitude based on default scaling relationship.

    Parameters
    ----------
    magnitude : float
        The magnitude of the rupture.

    Returns
    -------
    float
        An estimate of PGV for the rupture.
    """
    return np.interp(
        magnitude,
        [3.5, 4.1, 4.7, 5.2, 5.5, 5.8, 6.2, 6.5, 6.8, 7.0, 7.4, 7.7, 8.0],
        [0.015, 0.0375, 0.075, 0.15, 0.25, 0.4, 0.7, 1.0, 1.35, 1.65, 2.1, 2.5, 3.0],
    )


def find_rrup(magnitude: float, avg_dip: float, avg_rake: float) -> Tuple[float, float]:
    """
    Find rrup at which pgv estimated from magnitude is close to target.

    Parameters
    ----------
    pgv_target : pgv that we wish to come close to

    Returns
    -------
    Tuple[float, float]
        The (rrup, pgv) pair.
    """

    pgv_target = pgv_estimate_from_magnitude(magnitude)

    def pgv_delta_from_rrup(rrup: float):
        # stolen from rel2vm_params
        siteprop = Site()
        siteprop.vs30 = 500
        siteprop.Rtvz = 0
        siteprop.vs30measured = False
        siteprop.z1p0 = estimate_z1p0(siteprop.vs30)
        faultprop = Fault()
        faultprop.dip = avg_dip
        faultprop.rake = avg_rake
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
    if pgv_delta > 1e-4:
        raise ValueError(f"Failed to converge on rrup optimisation.")
    return rrup, pgv_target + pgv_delta


def minimum_area_bounding_box_for_polygons_masked(
    polygons: list[Polygon], mask: Polygon
) -> BoundingBox:
    polygon_union = shapely.normalize(shapely.union_all(polygons))
    bounding_polygon = shapely.normalize(shapely.intersection(polygon_union, mask))
    if isinstance(bounding_polygon, shapely.Polygon):
        return bounding_box.minimum_area_bounding_box(
            np.array(bounding_polygon.exterior.coords)
        )
    return bounding_box.minimum_area_bounding_box(
        np.vstack([np.array(geom.exterior.coords) for geom in bounding_polygon.geoms])
    )


def guess_simulation_duration(
    bounding_box: BoundingBox,
    type5_realisation: realisation.Realisation,
    ds_multiplier: float,
) -> float:
    """Guess the simulation duration for a realisation simulated in a given domain.

    The simulation duration is calculated as the time for the s-waves of
    a rupture to propogate from the centre of the domain to the edge of
    the domain.

    Parameters
    ----------
    bounding_box : BoundingBox
        The simulation domain.
    type5_realisation : realisation.Realisation
        The realisation.
    ds_multiplier : float
        A multiplier for the wavelength of the s-wave.

    Returns
    -------
    float
        An estimated simulation duration time.
    """
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


def get_max_depth(magnitude: float, hypocentre_depth: float):
    """
    Estimate the maximum depth to simulate for a rupture at a given depth
    with a given magnitude.

    Parameters
    ----------
    magnitude : float
        The magnitude of the rupture.
    hypocentre_depth : float
        hypocentre depth.


    Returns
    -------
    float
        The "interesting" simulation depth.
    """
    return round(
        10
        + hypocentre_depth
        + (
            10
            * np.power(
                (0.5 * np.power(10, (0.55 * magnitude - 1.2)) / hypocentre_depth), 0.3
            )
        ),
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
    min_vs: Annotated[
        float,
        typer.Option(
            help="Minimum velocity (in km/s),"
            "used to determine the boundary between LF and HF calculation."
        ),
    ] = 0.5,
    ds_multiplier: Annotated[float, typer.Option(help="Ds multiplier")] = 1.2,
    vm_version: Annotated[str, typer.Option(help="Velocity model version.")] = "2.06",
    vm_topo_type: Annotated[
        str, typer.Option(help="VM topology type")
    ] = "SQUASHED_TAPERED",
):
    "Generate velocity model parameters for a Type-5 realisation."
    type5_realisation = realisation.read_realisation(realisation_filepath)
    minimum_bounding_box = bounding_box.minimum_area_bounding_box(
        np.vstack(
            [fault.corners_nztm() for fault in type5_realisation.faults.values()]
        )[:, :2]
    )
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
    (rrup, _) = find_rrup(
        type5_realisation.magnitude,
        np.mean([fault.planes[0].dip for fault in type5_realisation.faults.values()]),
        np.mean([fault.planes[0].rake for fault in type5_realisation.faults.values()]),
    )
    site_inclusion_polygon = shapely.Point(minimum_bounding_box.origin).buffer(
        rrup * 1000
    )
    optimal_bounding_box = minimum_area_bounding_box_for_polygons_masked(
        [minimum_bounding_box.polygon, site_inclusion_polygon], get_land_polygon()
    )
    nx = int(np.ceil(optimal_bounding_box.extent_x / resolution))
    ny = int(np.ceil(optimal_bounding_box.extent_y / resolution))
    nz = int(np.ceil((max_depth - min_depth) / (resolution)))
    sim_duration = guess_simulation_duration(
        optimal_bounding_box, type5_realisation, ds_multiplier
    )
    box_origin_coords = coordinates.nztm_to_wgs_depth(
        np.append(optimal_bounding_box.origin, 0)
    )
    normalised_realisation = realisation.normalise_name(type5_realisation.name)
    vm_params = {
        "mag": magnitude,
        "MODEL_LAT": float(box_origin_coords[0]),
        "MODEL_LON": float(box_origin_coords[1]),
        "MODEL_ROT": float(optimal_bounding_box.bearing),
        "hh": resolution,
        "extent_x": float(optimal_bounding_box.extent_x),
        "extent_y": float(optimal_bounding_box.extent_y),
        "extent_zmax": float(max_depth),
        "extent_zmin": float(min_depth),
        "sim_duration": float(sim_duration),
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "min_vs": min_vs,
        "flo": min_vs / (5.0 * resolution),
        "sufx": f"_{normalised_realisation}01-h{resolution:.3f}",
        "model_version": vm_version,
        "topo_type": vm_topo_type,
    }
    with open(output_path, "w", encoding="utf-8") as output_file_handle:
        yaml.dump(vm_params, output_file_handle)


if __name__ == "__main__":
    typer.run(main)

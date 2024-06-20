#!/usr/bin/env python3
"""
VM Parameters Generation

This script generates the velocity model parameters used to generate the velocity model.

Usage
-----
To generate VM parameters for a Type-5 realisation:

```
$ python vm_params_generation.py path/to/realisation.yaml output/vm_params.yaml
```
"""


from pathlib import Path
from typing import Annotated, Tuple

import numpy as np
import pandas as pd
import scipy as sp
import shapely
import typer
import yaml
from empirical.util import openquake_wrapper_vectorized as openquake
from empirical.util import z_model_calculations
from empirical.util.classdef import GMM, TectType
from models import AfshariStewart_2016_Ds, classdef
from qcore import bounding_box, coordinates
from qcore.bounding_box import BoundingBox
from shapely import Polygon

from srf_generation import realisation

script_dir = Path(__file__).resolve().parent
NZ_LAND_OUTLINE = script_dir / "../SrfGen/NHM/res/rough_land.txt"


def get_nz_outline_polygon() -> Polygon:
    """
    Get the outline polygon of New Zealand.

    Returns:
        Polygon: The outline polygon of New Zealand.
    """
    polygon_coordinates = np.loadtxt(NZ_LAND_OUTLINE)
    polygon_coordinates = polygon_coordinates[:, ::-1]
    polygon_coordinates = np.append(
        polygon_coordinates,
        np.zeros_like(polygon_coordinates[:, 0]).reshape((-1, 1)),
        axis=1,
    )
    polygon_coordinates = coordinates.wgs_depth_to_nztm(polygon_coordinates)[:, :2]
    return shapely.Polygon(polygon_coordinates)


def pgv_estimate_from_magnitude(magnitude: np.ndarray) -> np.ndarray:
    """
    Return PGV for a given magnitude based on default scaling relationship.

    Parameters
    ----------
    magnitude : np.ndarray
        The magnitude(s) of the rupture(s).

    Returns
    -------
    np.ndarray
        An estimate of PGV for the rupture(s).

    References
    ----------
    See the "Custom Models Used in VM Params" wiki page for an explanation of this function.
    """
    return np.interp(
        magnitude,
        [3.5, 4.1, 4.7, 5.2, 5.5, 5.8, 6.2, 6.5, 6.8, 7.0, 7.4, 7.7, 8.0],
        [0.015, 0.0375, 0.075, 0.15, 0.25, 0.4, 0.7, 1.0, 1.35, 1.65, 2.1, 2.5, 3.0],
    )


def find_rrup(magnitude: float, avg_dip: float, avg_rake: float) -> Tuple[float, float]:
    """Find rrup at which pgv estimated from magnitude is close to target.

    Estimates rrup by calculating the rrup value that produces an
    estimated PGV value when measured from a site rrup distance
    away. The pgv -> rrup estimation comes from Chiou and Young (2014) [0].

    Parameters
    ----------
    magnitude : float
        The magnitude of the rupture.
    avg_dip : float
        The average dip of the faults involved.
    avg_rake : float
        The average rake of the faults involved.

    Returns
    -------
    Tuple[float, float]
        The (rrup, pgv) pair.

    References
    ----------
    [0]: Chiou BS-J, Youngs RR. Update of the
    Chiou and Youngs NGA Model for the Average Horizontal Component of
    Peak Ground Motion and Response Spectra. Earthquake
    Spectra. 2014;30(3):1117-1153.
    """

    pgv_target = pgv_estimate_from_magnitude(magnitude)

    def pgv_delta_from_rrup(rrup: float):
        vs30 = 500
        oq_dataframe = pd.DataFrame.from_dict(
            {
                "vs30": [vs30],
                "vs30measured": [False],
                "z1pt0": [z_model_calculations.chiou_young_08_calc_z1p0(vs30)],
                "dip": [avg_dip],
                "rake": [avg_rake],
                "mag": [magnitude],
                "ztor": [0],
                "rrup": [rrup],
                "rx": [rrup],
                "rjb": [rrup],
            }
        )
        pgv = openquake.oq_run(
            GMM.CY_14,
            TectType.ACTIVE_SHALLOW,
            oq_dataframe,
            "PGV",
        )["PGV_mean"].iloc[0]
        return np.abs(pgv - pgv_target)

    rrup_optimise_result = sp.optimize.minimize_scalar(
        pgv_delta_from_rrup, bounds=(0, 1e4)
    )
    rrup = rrup_optimise_result.x
    pgv_delta = rrup_optimise_result.fun
    if pgv_delta > 1e-4:
        raise ValueError("Failed to converge on rrup optimisation.")
    return rrup, pgv_target + pgv_delta


def estimate_simulation_duration(
    bounding_box: BoundingBox,
    type5_realisation: realisation.Realisation,
    ds_multiplier: float,
) -> float:
    """Estimate the simulation duration for a realisation simulated in a given domain.

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
    siteprop = classdef.Site()
    siteprop.vs30 = 500
    siteprop.Rtvz = 0
    siteprop.vs30measured = False
    siteprop.z1p0 = classdef.estimate_z1p0(siteprop.vs30)
    siteprop.Rrup = largest_corner_distance
    faultprop = classdef.Fault()
    faultprop.ztor = (
        0.0  # DON'T CHANGE THIS - this assumes we have a surface point source
    )
    faultprop.tect_type = classdef.TectType.ACTIVE_SHALLOW
    faultprop.faultstyle = classdef.FaultStyle.UNKNOWN
    faultprop.Mw = type5_realisation.magnitude

    ds = AfshariStewart_2016_Ds.Afshari_Stewart_2016_Ds(siteprop, faultprop, "Ds595")[0]

    return s_wave_arrival_time + ds_multiplier * ds


def get_max_depth(magnitude: float, hypocentre_depth: float) -> int:
    """
    Estimate the maximum depth to simulate for a rupture at a given depth
    with a given magnitude.

    Parameters
    ----------
    magnitude : float
        The magnitude of the rupture.
    hypocentre_depth : float
        hypocentre depth (in km).

    Returns
    -------
    float
        The maximum simulation depth.

    References
    ----------
    See the "Custom Models Used in VM Params" wiki page for an explanation of this function.
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

    initial_fault = type5_realisation.initial_fault()
    magnitude = type5_realisation.magnitude
    min_depth = 0

    # Get max depth and rrup
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

    # Get bounding box
    site_inclusion_polygon = shapely.Point(minimum_bounding_box.origin).buffer(
        rrup * 1000
    )
    optimal_bounding_box = bounding_box.minimum_area_bounding_box_for_polygons_masked(
        [minimum_bounding_box.polygon],
        [site_inclusion_polygon],
        get_nz_outline_polygon(),
    )

    # Calculate velocity model discretisation parameters
    nx = int(np.ceil(optimal_bounding_box.extent_x / resolution))
    ny = int(np.ceil(optimal_bounding_box.extent_y / resolution))
    nz = int(np.ceil((max_depth - min_depth) / (resolution)))
    sim_duration = estimate_simulation_duration(
        optimal_bounding_box, type5_realisation, ds_multiplier
    )

    box_origin_coords = coordinates.nztm_to_wgs_depth(
        np.append(optimal_bounding_box.origin, 0)
    )

    # Write the VM parameters file
    normalised_realisation_name = realisation.normalise_name(type5_realisation.name)

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
        "sufx": f"_{normalised_realisation_name}01-h{resolution:.3f}",
        "model_version": vm_version,
        "topo_type": vm_topo_type,
    }

    with open(output_path, "w", encoding="utf-8") as output_file_handle:
        yaml.dump(vm_params, output_file_handle)


if __name__ == "__main__":
    typer.run(main)

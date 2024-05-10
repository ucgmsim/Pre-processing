#!/usr/bin/env python3
"""
Generate GeoJSON objects for visualising type-5 rupture realisations.

This script converts rupture realisations into GeoJSON format for visualisation.

Usage:
    visualise_realisation.py REALISATION_FILE

Arguments:
    REALISATION_FILE  The YAML file containing a type-5 realisation.

Options:
    -h, --help  Show this help message and exit.

Examples:
    visualise_realisation.py my_realisation.yaml
"""

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from srf_generation import fault, realisation


def count_prior_triggered_faults(fault: fault.Fault) -> int:
    """
    Return the number of faults that trigger prior to the given fault.

    Parameters
    ----------
    fault : fault.Fault
        The fault object for which to determine the triggering level.

    Returns
    -------
    int
        The number of faults that trigger prior to the given fault.
    """
    i = 0
    cur = fault
    while cur.parent is not None:
        cur = cur.parent
        i += 1
    return i


def generate_darkening_colors(num_colors: int) -> list[(int, int, int)]:
    """
    Generate a list of darkening colors from red to black.

    Parameters
    ----------
    num_colors : int
        The number of colors to generate.

    Returns
    -------
    list[(int, int, int)]
        A list of RGB tuples representing darkening colors beginning with
        #FF0000 and ending at #000000 in linear steps.

    """
    colors = []
    r, g, b = 255, 0, 0

    step = 255 // num_colors

    for _ in range(num_colors):
        colors.append((r, g, b))
        r = max(r - step, 0)
        g = max(g - step, 0)
        b = max(b - step, 0)

    return colors


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """
    Convert RGB tuple to hex string.

    Parameters
    ----------
    r : int
        The red channel value (0 to 255).
    g : int
        The green channel value (0 to 255).
    b : int
        The blue channel value (0 to 255).

    Returns
    -------
    str
        The hex string representation of the RGB color.

    """
    return "#{:02x}{:02x}{:02x}".format(r, g, b)


def fault_to_geojson_object(
    fault: fault.Fault, fault_colour: (int, int, int)
) -> dict[str, Any]:
    """
    Convert Fault object to GeoJSON object.

    Parameters
    ----------
    fault : fault.Fault
        The Fault object to convert.
    fault_colour : (int, int, int)
        The RGB tuple representing the color of the fault.

    Returns
    -------
    dict
        The GeoJSON representation of the fault geometry.

    """
    geojson_object = {
        "type": "Feature",
        "properties": {
            "fill": rgb_to_hex(*fault_colour),
            "fill-opacity": 1,
            "name": fault.name,
            "triggered by": (
                "initial fault" if fault.parent is None else fault.parent.name
            ),
            "magnitude": fault.magnitude,
        },
    }
    if fault.segments[0].dip == 90:
        corners = fault.corners()[:, :2, :2].reshape((-1, 2))[:, ::-1]
        geojson_object["geometry"] = {
            "type": "LineString",
            "coordinates": corners.tolist(),
        }
    else:
        corners = fault.corners()[:, :, :2]
        first_corner = corners[:, 0, :].reshape((-1, 1, 2))
        polygon_corners = np.concatenate((corners, first_corner), axis=1).reshape(
            -1, 1, 5, 2
        )[:, :, :, ::-1]
        geojson_object["geometry"] = {
            "type": "MultiPolygon",
            "coordinates": polygon_corners.tolist(),
        }

    return geojson_object


def fault_hypocentre_geojson(fault: fault.Fault) -> dict[str, Any]:
    """
    Generate GeoJSON representation of fault hypocentre.

    Parameters
    ----------
    fault : fault.Fault
        The Fault object containing hypocentre information.

    Returns
    -------
    dict[str, Any]
        The GeoJSON representation of the fault hypocentre.

    """
    geojson_object = {
        "type": "Feature",
        "geometry": {
            "type": "Point",
            "coordinates": fault.fault_coordinates_to_wgsdepth_coordinates(
                np.array([fault.shyp, fault.dhyp])
            )[:2][::-1].tolist(),
        },
        "properties": {"name": "hypocentre", "marker-color": "red"},
    }
    return geojson_object


def fault_jump_point_geojson(fault: fault.Fault) -> dict[str, Any]:
    """
    Generate GeoJSON representation of fault jump point.

    Parameters
    ----------
    fault : fault.Fault
        The Fault object containing jump point information.

    Returns
    -------
    dict[str, Any]
        The GeoJSON representation of the fault jump point.

    """
    geojson_object = {
        "type": "Feature",
        "geometry": {
            "type": "Point",
            "coordinates": fault.parent_jump_coords[:2][::-1],
        },
        "properties": {"name": "jump point", "marker-colour": "blue"},
    }
    return geojson_object


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="visualise_realisation",
        description="Generate a geojson object to visualise a type-5 realisation (for sanity checking).",
    )
    parser.add_argument(
        "realisation_file",
        type=Path,
        metavar="REALISATION_FILE",
        help="The YAML file containing a type-5 realisation",
    )
    args = parser.parse_args()
    realisation = realisation.read_realisation(args.realisation_file)

    max_fault_level = (
        max(
            count_prior_triggered_faults(fault) for fault in realisation.faults.values()
        )
        + 1
    )

    colours = generate_darkening_colors(max_fault_level)
    features = []
    geojson_object = {"type": "FeatureCollection", "features": features}

    for fault in realisation.faults.values():
        level = count_prior_triggered_faults(fault)
        features.append(fault_to_geojson_object(fault, colours[level]))
        if fault.shyp is not None and fault.dhyp is not None:
            features.append(fault_hypocentre_geojson(fault))
            if fault.parent:
                features.append(fault_jump_point_geojson(fault))

    print(json.dumps(geojson_object))

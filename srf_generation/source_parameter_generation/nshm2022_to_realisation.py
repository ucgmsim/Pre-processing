#!/usr/bin/env python3
"""
Generate realisation stub files from ruptures in the NSHM 2022 database.

This script generates YAML realisation stub files from ruptures in the NSHM 2022
database. It extracts fault geometry computes causality information from the database
and incorporates default parameter values to generate the realisation.

Usage:
    nshm2022_to_realisation.py NSHM_DB_FILE RUPTURE_ID YAML_FILE [--dt DT] [--genslip-seed GENSLIP_SEED] [--srfgen-seed SRFGEN_SEED]

Arguments:
    RUPTURE_ID      The ID of the rupture to generate the realisation stub for
                    (find this using the NSHM Rupture Explorer).
    YAML_FILE       Location to write out the YAML realisation value.
    NSHM_DB_FILE    The NSHM sqlite database containing rupture
                    information and fault geometry.

Options:
    --help                      Show help message and exit.
    --dt DT                     Time resolution for source modelling (default:
                                0.0500s).
    --genslip-seed GENSLIP_SEED Seed for genslip, used to initialise slip
                                distribution on fault (default: 1).
    --srfgen-seed SRFGEN_SEED   Seed for srfgen, used to initialise slip
                                distribution on fault (default: 1).

Depends On:
- python >= 3.11
- numpy
- pyproj
- qcore
- scipy
- typer
- PyYAML
"""
import argparse
import collections
import itertools
import re
import sqlite3
from pathlib import Path
from sqlite3 import Connection
from typing import Annotated, Any, TextIO

import numpy as np
import pyproj
import qcore.geo
import qcore.uncertainties.mag_scaling
import scipy as sp
import typer
import yaml
from srf_generation import fault

import rupture_propogation
from rupture_propogation import RuptureCausalityTree

WGS_CODE = 4326
NZTM_CODE = 2193
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)

app = typer.Typer()


def get_faults_for_rupture(
    connection: sqlite3.Connection, rupture_id: int
) -> list[fault.Fault]:
    """Get all fault sections that participate in a given rupture.

    Parameters
    ----------
    connection : sqlite3.Connection
        The SQLite connection of the NSHM DB.
    rupture_id : int
        The rupture to query for.

    Returns
    -------
    list[fault.Fault]
        The list of faults involved in the rupture.
    """
    cursor = connection.cursor()
    cursor.execute(
        """SELECT fs.*, p.parent_id, p.name
        FROM fault_segment fs
        JOIN rupture_faults rf ON fs.fault_id = rf.fault_id
        JOIN fault f ON fs.fault_id = f.fault_id
        JOIN parent_fault p ON f.parent_id = p.parent_id
        WHERE rf.rupture_id = ? AND fs.length >= 0.1
        ORDER BY f.parent_id""",
        (rupture_id,),
    )

    fault_segments = cursor.fetchall()
    cur_parent_id = None
    faults = []
    for (
        _,
        strike,
        rake,
        dip,
        dtop,
        dbottom,
        length,
        width,
        dip_dir,
        clon,
        clat,
        _,
        parent_id,
        parent_name,
    ) in fault_segments:
        if parent_id != cur_parent_id:
            # NOTE: The location of the hypocentre for this fault in fault
            # coordinates (that is: shyp and dhyp) are determined dynamically
            # based on the rupture causality tree.
            faults.append(
                fault.Fault(
                    name=parent_name, tect_type=None, segments=[], dhyp=None, shyp=None
                )
            )
            cur_parent_id = parent_id
        segment = fault.FaultSegment(
            strike, rake, dip, dtop, dbottom, length, width, dip_dir, clon, clat
        )
        faults[-1].segments.append(segment)
    return faults


def wgsdepth_to_nztm(wgsdepthcoordinates: np.ndarray) -> np.ndarray:
    """Convert (lat, lon, depth) coordinates to NZTM coordinates.

    Parameters
    ----------
    wgsdepthcoordinates : np.ndarray
        An (n x 3) array of (lat, lon, depth) coordinates.

    Returns
    -------
    np.ndarray
        An (n x 3) array of (x, y, depth) coordinates, where x and y specify
        coordinates in the NZTM coordinate system.
    """
    nztm_coords = np.array(
        WGS2NZTM.transform(wgsdepthcoordinates[:, 0], wgsdepthcoordinates[:, 1]),
    ).T
    return np.append(nztm_coords, wgsdepthcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def nztm_to_wgsdepth(nztmcoordinates: np.ndarray) -> np.ndarray:
    """Convert NZTM coordinates to WGS84+depth coordinates

    Parameters
    ----------
    nztmcoordinates : np.ndarray
        An (n x 3) array of NZTM coordinates (x, y, depth).

    Returns
    -------
    np.ndarray
        An (n x 3) array of WGS84 coordinates (lat, lon, depth).
    """
    wgs_coords = np.array(
        NZTM2WGS.transform(nztmcoordinates[:, 0], nztmcoordinates[:, 1]),
    ).T
    return np.append(wgs_coords, nztmcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def closest_points_between_faults(
    fault_u: fault.Fault, fault_v: fault.Fault
) -> (np.ndarray, np.ndarray):
    """Given two faults u and v, return the closest pair of points on the faults.

    Parameters
    ----------
    fault_u : fault.Fault
        The fault u.
    fault_v : fault.Fault
        The fault v.

    Returns
    -------
    np.ndarray
        The pair of points (x, y) with x in fault u and y in fault such that the
        distance (in 3 dimensions) between x and y is minimised. The points x
        and y are given in the WGS84 coordinate system.
    """
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
    """Given two fault segments, establish if the faults could possibly be closer than cutoff distance from each other.

    This function is used to filter out far away fault segments, because
    computing the exact closest distance is slow.

    Parameters
    ----------
    fault1 : fault.FaultSegment
        The fault u.
    fault2 : fault.FaultSegment
        The fault v.
    cutoff : float
        The cutoff distance (in metres).

    Returns
    -------
    bool
        Returns True if it is possible that the fault segments fault1 and fault2
        could be within cutoff metres from each other at their closest distance.
    """
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
    """Given two faults, establish if any fault segments could be closer than cutoff distance from each other.

    This function is used to filter out far apart faults, because computing the
    exact closest distance is slow.

    Parameters
    ----------
    fault1 : fault.Fault
        The fault u.
    fault2 : fault.Fault
        The fault v.
    cutoff : float
        The cutoff distance (in metres).

    Returns
    -------
    bool
        Return True if any pair of fault segments from u and v could be closer
        than cutoff metres from each other.
    """
    return any(
        test_fault_segment_vialibility(seg1, seg2, cutoff)
        for (seg1, seg2) in itertools.product(fault1.segments, fault2.segments)
    )


def build_rupture_causality_tree(
    initial_fault: fault.Fault, faults: list[fault.Fault]
) -> RuptureCausalityTree:
    """Compute the rupture causality tree for a given series of faults, rupturing from a given initial fault.

    Parameters
    ----------
    initial_fault : fault.Fault
        The initial fault to rupture from.
    faults : list[fault.Fault]
        The list of faults that must rupture, possible including the initial
        fault.

    Returns
    -------
    RuptureCausalityTree
        A tree mapping faults to their parent (triggering) fault. The initial
        fault has the parent None.
    """
    jump_point_map = collections.defaultdict(dict)
    cutoff = 15000
    for fault_u in faults:
        fault_u_name = fault_u.name
        for fault_v in faults:
            fault_v_name = fault_v.name
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
    cutoff = 15000
    pruned = rupture_propogation.prune_distance_graph(distance_graph, cutoff)
    probability_graph = rupture_propogation.probability_graph(pruned)

    return rupture_propogation.probabilistic_minimum_spanning_tree(
        probability_graph, initial_fault.name
    )


def compute_jump_point_hypocentre_fault_coordinates(
    from_fault: fault.Fault, to_fault: fault.Fault
) -> (np.ndarray, np.ndarray):
    """Compute the closest two points between two faults.

    Parameters
    ----------
    from_fault : fault.Fault
        The fault to jump from.
    to_fault : fault.Fault
        The fault to jump to.

    Returns
    -------
    (np.ndarray, np.ndarray)
        The closest pair of points (u, v) with u in from_fault and v in
        to_fault. The points are specified in WGS84+depth coordinates.
    """

    from_fault_point_nztm, to_fault_point_nztm = closest_points_between_faults(
        from_fault, to_fault
    )
    from_fault_point_wgsdepth = nztm_to_wgsdepth(
        from_fault_point_nztm.reshape((1, -1))
    ).ravel()
    to_fault_point_wgsdepth = nztm_to_wgsdepth(
        to_fault_point_nztm.reshape((1, -1))
    ).ravel()
    return from_fault_point_wgsdepth, to_fault_point_wgsdepth


def link_hypocentres(
    rupture_causality_tree: RuptureCausalityTree, faults: list[fault.Fault]
):
    """Set the jumping points across faults.

    Given a rupture causality tree, and the list of faults on that tree, set the
    parent, shyp, dhyp and jump point coordinates based on the computed fault
    jumping points. Modifies the faults array.

    Parameters
    ----------
    rupture_causality_tree : RuptureCausalityTree
        The rupture causality tree that specifies how the faults jump.
    faults : list[fault.Fault]
        The list of faults to link.
    """
    fault_name_map = {fault.name: fault for fault in faults}
    for to_fault in faults:
        fault_name = to_fault.name
        if rupture_causality_tree[fault_name] is None:
            shyp, dhyp = to_fault.expected_fault_coordinates()
            to_fault.shyp = shyp
            to_fault.dhyp = dhyp
        else:
            from_fault = fault_name_map[rupture_causality_tree[fault_name]]
            from_fault_point, to_fault_point = (
                compute_jump_point_hypocentre_fault_coordinates(from_fault, to_fault)
            )
            to_shyp, to_dhyp = to_fault.hypocentre_wgs_to_fault_coordinates(
                to_fault_point
            )
            to_fault.parent = from_fault
            to_fault.shyp = float(to_shyp)
            to_fault.dhyp = float(to_dhyp)
            to_fault.parent_jump_coords = tuple(float(x) for x in from_fault_point)


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


def magnitude_for_fault(target_fault: fault.Fault, total_area: float) -> float:
    """Calculate the target fault's contribution to the total magnitude of a given rupture.

    Given a target fault, and the total area of a rupture the fault participates
    in, calculation the proportion of the magnitude that the target fault would
    produce when it ruptures.

    Parameters
    ----------
    target_fault : fault.Fault
        The fault to compute the target magnitude for.
    total_area : float
        The total area of the rupture the target fault participates in.

    Returns
    -------
    float
        The magnitude the fault must rupture at to participate.
    """
    log_area = np.log10(total_area)
    total_magnitude = qcore.uncertainties.mag_scaling.a_to_mw_leonard(
        total_area, 4.00, 3.99, 0
    )
    return float(np.log10(target_fault.area()) + total_magnitude - log_area)


def set_magnitudes(faults: fault.Fault):
    """Set the magnitude for each participating fault. Modifies the faults array.

    Parameters
    ----------
    faults : fault.Fault
        The faults participating in the rupture.
    """
    total_area = sum(fault.area() for fault in faults)
    for fault in faults:
        fault.magnitude = magnitude_for_fault(fault, total_area)


def write_yaml_realisation_stub_file(
    yaml_realisation_file: TextIO,
    default_parameter_values: dict[str, Any],
    faults: list[fault.Fault],
):
    """Write yaml realisation stub file.

    This function generates a yaml file stub for realisation, incorporating
    default parameter values, fault geometry, and rupture causality.

    Parameters
    ----------
    yaml_realisation_file : TextIO
        The file object to write the yaml to.
    default_parameter_values : dict[str, Any]
        Default parameter values for the realisation.
    faults : list[fault.Fault]
        The list of faults contaning both geoemetry and causality information.
    """
    realisation_object = {}
    realisation_object |= default_parameter_values
    fault_object = {}
    for fault in faults:
        segments = [
            {
                "strike": segment.strike,
                "rake": segment.rake,
                "dip": segment.dip,
                "dtop": segment.dtop,
                "dbottom": segment.dbottom,
                "length": segment.length,
                "width": segment.width,
                "dip_dir": segment.dip_dir,
                "clon": segment.clon,
                "clat": segment.clat,
            }
            for segment in fault.segments
        ]
        fault_object[normalise_name(fault.name)] = {
            "name": fault.name,
            "tect_type": fault.tect_type,
            "parent": normalise_name(fault.parent.name) if fault.parent else None,
            "magnitude": fault.magnitude,
            "shyp": fault.shyp,
            "dhyp": fault.dhyp,
            "parent_jump_coords": (
                list(fault.parent_jump_coords) if fault.parent else None
            ),
            "segments": segments,
        }
    realisation_object["faults"] = fault_object
    yaml.dump(realisation_object, yaml_realisation_file)


@app.command(
    help="Generate realisation stub files from ruptures in the NSHM 2022 database."
)
def main(
    nshm_db_file: Annotated[
        Path,
        typer.Argument(
            help="The NSHM sqlite database containing rupture information and fault geometry.",
            readable=True,
            exists=True,
        ),
    ],
    rupture_id: Annotated[
        int,
        typer.Argument(
            help="The ID of the rupture to generate the realisation stub for (find this using the NSHM Rupture Explorer)."
        ),
    ],
    yaml_file: Annotated[
        Path,
        typer.Argument(
            help="Location to write out the YAML realisation value.", writable=True
        ),
    ],
    dt: Annotated[
        float, typer.Option(help="Time resolution for source modelling.", min=0)
    ] = 0.05,
    genslip_seed: Annotated[
        int,
        typer.Option(
            help="Seed for genslip, used to initialise slip distribution on fault."
        ),
    ] = 1,
    srfgen_seed: Annotated[
        int,
        typer.Option(
            help="Seed for srfgen, used to initialise slip distribution on fault."
        ),
    ] = 1,
):
    with sqlite3.connect(nshm_db_file) as conn:
        faults = get_faults_for_rupture(conn, rupture_id)
    set_magnitudes(faults)
    initial_fault = faults[0]
    rupture_causality_tree = build_rupture_causality_tree(initial_fault, faults)
    link_hypocentres(rupture_causality_tree, faults)
    yaml_realisation_file = yaml_file
    default_parameter_values_with_args = {
        "name": f"Rupture {rupture_id}",
        "type": 5,
        "genslip_version": "5.4.2",
        "srfgen_seed": srfgen_seed,
        "genslip_seed": genslip_seed,
        "dt": dt,
    }
    with open(yaml_realisation_file, "w", encoding="utf-8") as yaml_out:
        write_yaml_realisation_stub_file(
            yaml_out, default_parameter_values_with_args, faults
        )


if __name__ == "__main__":
    app()

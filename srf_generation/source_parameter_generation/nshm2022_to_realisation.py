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
- nshmdb
- pyproj
- qcore
- scipy
- typer
- PyYAML
"""
import collections
import itertools
import re
from pathlib import Path
from typing import Annotated, Any, TextIO

import numpy as np
import qcore.coordinates
import qcore.geo
import qcore.uncertainties.mag_scaling
import rupture_propogation
import scipy as sp
import typer
import yaml
from nshmdb import nshmdb
from nshmdb.fault import FaultPlane
from rupture_propogation import RuptureCausalityTree
from srf_generation.realisation import RealisationFault

app = typer.Typer()


def closest_points_between_faults(
    fault_u: RealisationFault, fault_v: RealisationFault
) -> (np.ndarray, np.ndarray):
    """Given two faults u and v, return the closest pair of points on the faults.

    Parameters
    ----------
    fault_u : RealisationFault
        The fault u.
    fault_v : RealisationFault
        The fault v.

    Returns
    -------
    np.ndarray
        The pair of points (x, y) with x in fault u and y in fault such that the
        distance (in 3 dimensions) between x and y is minimised. The points x
        and y are given in the WGS84 coordinate system.
    """
    fault_u_corners = qcore.coordinates.wgs_depth_to_nztm(fault_u.corners()).reshape(
        (-1, 4, 3)
    )

    fault_v_corners = qcore.coordinates.wgs_depth_to_nztm(fault_v.corners()).reshape(
        (-1, 4, 3)
    )
    point_u, point_v = qcore.geo.closest_points_between_plane_sequences(
        fault_u_corners, fault_v_corners
    )
    return point_u, point_v


def test_fault_plane_vialibility(
    fault1: FaultPlane, fault2: FaultPlane, cutoff: float
) -> bool:
    """Given two fault planes, establish if the faults could possibly be closer than cutoff distance from each other.

    This function is used to filter out far away fault planes, because
    computing the exact closest distance is slow.

    Parameters
    ----------
    fault1 : FaultPlane
        The fault u.
    fault2 : FaultPlane
        The fault v.
    cutoff : float
        The cutoff distance (in metres).

    Returns
    -------
    bool
        Returns True if it is possible that the fault planes fault1 and fault2
        could be within cutoff metres from each other at their closest distance.
    """
    fault1_centroid = qcore.coordinates.wgs_depth_to_nztm(
        fault1.centroid().reshape((1, -1))
    ).ravel()
    fault2_centroid = qcore.coordinates.wgs_depth_to_nztm(
        fault2.centroid().reshape((1, -1))
    ).ravel()
    fault1_radius = max(fault1.width_m, fault1.length_m) + cutoff
    fault2_radius = max(fault2.width_m, fault2.length_m) + cutoff
    return qcore.geo.spheres_intersect(
        fault1_centroid, fault1_radius, fault2_centroid, fault2_radius
    )


def test_fault_viability(
    fault1: RealisationFault, fault2: RealisationFault, cutoff: float
) -> bool:
    """Given two faults, establish if any fault planes could be closer than cutoff distance from each other.

    This function is used to filter out far apart faults, because computing the
    exact closest distance is slow.

    Parameters
    ----------
    fault1 : Fault
        The fault u.
    fault2 : Fault
        The fault v.
    cutoff : float
        The cutoff distance (in metres).

    Returns
    -------
    bool
        Return True if any pair of fault planes from u and v could be closer
        than cutoff metres from each other.
    """
    return any(
        test_fault_plane_vialibility(seg1, seg2, cutoff)
        for (seg1, seg2) in itertools.product(fault1.planes, fault2.planes)
    )


def build_rupture_causality_tree(
    initial_fault: RealisationFault, faults: list[RealisationFault]
) -> RuptureCausalityTree:
    """Compute the rupture causality tree for a given series of faults, rupturing from a given initial fault.

    Parameters
    ----------
    initial_fault : RealisationFault
        The initial fault to rupture from.
    faults : list[RealisationFault]
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
    from_fault: RealisationFault, to_fault: RealisationFault
) -> (np.ndarray, np.ndarray):
    """Compute the closest two points between two faults.

    Parameters
    ----------
    from_fault : RealisationFault
        The fault to jump from.
    to_fault : RealisationFault
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
    from_fault_point_wgsdepth = qcore.coordinates.nztm_to_wgs_depth(
        from_fault_point_nztm.reshape((1, -1))
    ).ravel()
    to_fault_point_wgsdepth = qcore.coordinates.nztm_to_wgs_depth(
        to_fault_point_nztm.reshape((1, -1))
    ).ravel()
    return from_fault_point_wgsdepth, to_fault_point_wgsdepth


def link_hypocentres(
    rupture_causality_tree: RuptureCausalityTree, faults: list[RealisationFault]
):
    """Set the jumping points across faults.

    Given a rupture causality tree, and the list of faults on that tree, set the
    parent, shyp, dhyp and jump point coordinates based on the computed fault
    jumping points. Modifies the faults array.

    Parameters
    ----------
    rupture_causality_tree : RuptureCausalityTree
        The rupture causality tree that specifies how the faults jump.
    faults : list[RealisationFault]
        The list of faults to link.
    """
    fault_name_map = {fault.name: fault for fault in faults}
    for to_fault in faults:
        fault_name = to_fault.name
        if rupture_causality_tree[fault_name] is None:
            continue
        else:
            from_fault = fault_name_map[rupture_causality_tree[fault_name]]
            from_fault_point, to_fault_point = (
                compute_jump_point_hypocentre_fault_coordinates(from_fault, to_fault)
            )
            try:
                to_shyp, to_dhyp = to_fault.global_coordinates_to_fault_coordinates(
                    to_fault_point
                )
            except ValueError:
                breakpoint()
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


def magnitude_for_fault(target_fault: RealisationFault, total_area: float) -> float:
    """Calculate the target fault's contribution to the total magnitude of a given rupture.

    Given a target fault, and the total area of a rupture the fault participates
    in, calculation the proportion of the magnitude that the target fault would
    produce when it ruptures.

    Parameters
    ----------
    target_fault : RealisationFault
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


def set_magnitudes(faults: list[RealisationFault]):
    """Set the magnitude for each participating fault. Modifies the faults array.

    Parameters
    ----------
    faults : list[RealisationFault]
        The faults participating in the rupture.
    """
    total_area = sum(realisation_fault.area() for realisation_fault in faults)
    for realisation_fault in faults:
        realisation_fault.magnitude = magnitude_for_fault(realisation_fault, total_area)


def estimate_log_rupture_rate_for_magnitude(
    rupture_mfds: np.ndarray, magnitude: float
) -> float:

    def log_model(m, a, b):
        return a - b * m

    popt, _ = sp.optimize.curve_fit(
        log_model, rupture_mfds[:, 0], np.log(rupture_mfds[:, 1])
    )
    return log_model(magnitude, *popt)


def write_yaml_realisation_stub_file(
    yaml_realisation_file: TextIO,
    default_parameter_values: dict[str, Any],
    faults: list[RealisationFault],
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
    faults : list[RealisationFault]
        The list of faults contaning both geoemetry and causality information.
    """
    realisation_object = {}
    realisation_object |= default_parameter_values
    fault_object = {}
    for fault in faults:
        planes = [
            {"corners": plane.corners.tolist(), "rake": plane.rake}
            for plane in fault.planes
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
            "planes": planes,
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
    db = nshmdb.NSHMDB(nshm_db_file)
    faults = [
        RealisationFault(
            name=fault.name, tect_type=fault.tect_type, planes=fault.planes
        )
        for fault in db.get_rupture_faults(rupture_id)
    ]
    set_magnitudes(faults)
    initial_fault = faults[0]
    rupture_causality_tree = build_rupture_causality_tree(initial_fault, faults)
    link_hypocentres(rupture_causality_tree, faults)
    default_parameter_values_with_args = {
        "name": None,
        "type": 5,
        "genslip_version": "5.4.2",
        "srfgen_seed": srfgen_seed,
        "genslip_seed": genslip_seed,
        "dt": dt,
    }
    with open(yaml_file, "w", encoding="utf-8") as yaml_out:
        write_yaml_realisation_stub_file(
            yaml_out, default_parameter_values_with_args, faults
        )


if __name__ == "__main__":
    app()

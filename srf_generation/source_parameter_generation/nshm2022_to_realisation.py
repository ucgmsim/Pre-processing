#!/usr/bin/env python3
import sqlite3
import fault
import yaml
import rupture_prop

import rupture_propogation
from rupture_propogation import RuptureCausalityTree

def get_faults_for_rupture(connection: sqlite3.Connection, parent_map: dict[int, str], rupture_id: int) -> list[fault.Fault]
    cursor = connection.cursor()
    cursor.execute("""
        SELECT fs.*, f.parent_id, f.parent_name
        FROM fault_segment fs
        JOIN rupture_faults rf ON fs.segment_id = rf.segment_id
        JOIN fault f ON fs.fault_id = f.fault_id
        WHERE rf.rupture_id = ?
        ORDER BY f.parent_id
    """, (rupture_id,))

    fault_segments = cursor.fetchall()
    cur_parent_id = None
    faults = []
    for _, strike, rake, dip, dtop, dbottom, length, width, dip_dir, clon, clat, _, parent_id, parent_name in fault_segments:
        if parent_id != cur_parent_id:
            faults.append(Fault(name=parent_name, tect_type=None, segments=[], dhyp=None, shyp=None))
            cur_parent_id = parent_id
        segment = FaultSegment(strike, rake, dip, dtop, dbottom, length, width, dip_dir, clon, clat)
        faults[-1].segments.append(segment)
    return faults

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


def build_rupture_causality_tree(initial_fault: fault.Fault, faults: list[fault.Fault]) -> RuptureCausalityTree:
    jump_point_map = collections.defaultdict(dict) .
    cutoff = 15000
    for fault_u in realisation.faults.items():
        fault_u_name = fault_u.name
        for fault_v in realisation.faults.items():
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

    pruned = rupture_propogation.prune_distance_graph(distance_graph, 15000)
    probability_graph = rupture_propogation.probability_graph(pruned)

    return rupture_propogation.probabilistic_minimum_spanning_tree(
        probability_graph, initial_fault.name
    )


def compute_jump_point_hypocentre_fault_coordinates(
    from_fault: fault.Fault, to_fault: fault.Fault
) -> (float, float):
    from_fault_point_nztm, to_fault_point_nztm = closest_points_between_faults(from_fault, to_fault)
    to_fault_point_wgsdepth = to_fault.hypocentre_wgs_to_fault_coordinates(
        nztm_to_wgsdepth(to_fault_point_nztm.reshape((1, -1))).ravel()
    )
    from_fault_point_wgsdepth = from_fault.hypocentre_wgs_to_fault_coordinates(
        nztm_to_wgsdepth(from_fault_point_nztm.reshape((1, -1))).ravel()
    )
    return from_fault, to_fault


def link_hypocentres(
    rupture_causality_tree: RuptureCausalityTree
, faults: list[fault.Fault]):
    fault_name_map = {fault.name: fault for fault in faults}
    for to_fault in faults:
        fault_name = to_fault.name
        if rupture_causality_tree[fault_name] is None:
            continue
        else:
            from_fault = fault_name_map[rupture_causality_tree[fault_name]]
            (from_jump_strike, from_jump_dip), ( to_shyp, to_dhyp ) = compute_jump_point_hypocentre_fault_coordinates(from_fault, to_fault)
            to_fault.parent = from_fault
            to_fault.shyp = shyp
            to_fault.dhyp = dhyp
            to_fault.parent_jump_strike = from_jump_strike
            to_fault.parent_jump_dip = from_jump_dip


def write_yaml_realisation_stub_file(yaml_realisation_file, default_parameter_values, rupture_causality_tree, faults):
    realisation_object = {}
    realisation_object |= default_parameter_values
    fault_object = []
    for fault in faults:
        segments = [
            {
                'strike': segment.strike,
                'rake': segment.rake,
                'dip': segment.dip,
                'dtop': segment.dtop,
                'dbottom': segment.dbottom,
                'length': segment.length,
                'width': segment.width,
                'dip_dir': segment.dip_dir,
                'clon': segment.clon,
                'clat': segment.clat
            }
        ]
        fault_object.append({
            'name': fault.name,
            'tect_type': fault.tect_type
            'parent': fault.parent.name,
            'shyp': fault.shyp,
            'dhyp': fault.dhyp
            'parent_jump_strike': fault.parent_jump_strike,
            'parent_jump_dip': fault.parent_jump_dip,
            'segments': segments
        })
    yaml.dump(yaml_realisation_file, fault_object)




if __name__ == "__main__":
    rupture_id = 100
    with sqlite3.connect(nshm_db_file) as conn:
        faults = get_parent_faults_for_rupture(conn, rupture_id)
    initial_fault = faults[0]
    rupture_causality_tree = build_rupture_causality_tree(initial_fault, faults)
    link_hypocentres(rupture_causality_tree, faults)
    with open(yaml_realisation_file, 'w') as f:
        write_yaml_realisation_stub_file(yaml_realisation_file, default_parameter_values, faults)

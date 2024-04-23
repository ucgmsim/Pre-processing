#!/usr/bin/env python3
import json
import yaml
import fault
import qcore.geo
import numpy as np
import collections
import re


def strike_between_coordinates(a, b) -> float:
    a_lat, a_lon = a
    b_lat, b_lon = b
    return qcore.geo.ll_bearing(a_lon, a_lat, b_lon, b_lat)


def distance_between(a, b) -> float:
    a_lat, a_lon = a
    b_lat, b_lon = b
    return qcore.geo.ll_dist(a_lon, a_lat, b_lon, b_lat)


def centre_point(a, b, dip, dip_dir, width):
    a_lat, a_lon = a
    b_lat, b_lon = b
    c_lon, c_lat = qcore.geo.ll_mid(a_lon, a_lat, b_lon, b_lat)
    projected_width = width * np.cos(np.radians(dip)) / 2
    return qcore.geo.ll_shift(c_lat, c_lon, projected_width, dip_dir)


def geojson_feature_to_fault(feature):
    properties = feature["properties"]
    dip = properties["DipDeg"]
    rake = properties["Rake"]
    dbottom = properties["LowDepth"]
    dtop = properties["UpDepth"]
    dip_dir = properties["DipDir"]
    name = properties["ParentName"]
    leading_edge = feature["geometry"]
    width = float(dbottom / np.sin(np.radians(dip)))
    segments = []
    for i in range(len(leading_edge) - 1):
        left = tuple(reversed(leading_edge[i]))
        right = tuple(reversed(leading_edge[i + 1]))
        c_lat, c_lon = centre_point(left, right, dip, dip_dir, dbottom)
        strike = strike_between_coordinates(left, right)
        length = distance_between(left, right)
        segments.append(
            fault.FaultSegment(
                strike=strike,
                rake=rake,
                dip=dip,
                dtop=dtop,
                dbottom=dbottom,
                length=length,
                width=width,
                dip_dir=dip_dir,
                clon=c_lon,
                clat=c_lat,
            )
        )
    return fault.Fault(name=name, tect_type=None, segments=segments)


def merge_by_parents(fault_features):
    merged = {}
    for feature in fault_features:
        parent_name = feature["properties"]["ParentName"]
        if parent_name not in merged:
            merged[parent_name] = {
                "properties": feature["properties"],
                "geometry": feature["geometry"]["coordinates"],
            }
        else:
            merged[parent_name]["geometry"].extend(
                feature["geometry"]["coordinates"][1:]
            )
    return merged.values()


def normalise_name(name: str) -> str:
    fault_no_illegal_characters = re.sub(r" |,|:", "_", name)
    return fault_no_illegal_characters.lower()


CAUSALITY_MAP = {
    "alpine__caswell": None,
    "alpine__caswell_-_south_george": "alpine__caswell",
    "alpine__george_landward": "alpine__caswell_-_south_george",
    "alpine__george_to_jacksons": "alpine__george_landward",
    "alpine__jacksons_to_kaniere": "alpine__george_to_jacksons",
    "alpine__kaniere_to_springs_junction": "alpine__jacksons_to_kaniere",
    "alpine__nancy": "alpine__caswell",
    "alpine__resolution_-_dagg": "alpine__nancy",
    "alpine__resolution_-_five_fingers": "alpine__resolution_-_dagg",
    "booboo": "needles",
    "campbell_bank": "chancet",
    "chancet": "booboo",
    "clarence__central": "clarence__southwest",
    "clarence__southwest": "alpine__kaniere_to_springs_junction",
    "fidget": "clarence__central",
    "kekerengu_1": "fidget",
    "needles": "kekerengu_1",
}


def tree_level(tree, node):
    i = 0
    cur = node
    while tree[cur] is not None:
        cur = tree[cur]
        i += 1
    return i


def generate_darkening_colors(num_colors):
    colors = []
    # Start with fully saturated red (255, 0, 0)
    r, g, b = 255, 0, 0

    # Calculate the step size for each color channel
    step = 255 // num_colors

    for _ in range(num_colors):
        colors.append((r, g, b))
        # Decrease the intensity of each color channel
        r = max(r - step, 0)
        g = max(g - step, 0)
        b = max(b - step, 0)

    return colors


def rgb_to_hex(r, g, b):
    """Convert RGB tuple to hex string."""
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

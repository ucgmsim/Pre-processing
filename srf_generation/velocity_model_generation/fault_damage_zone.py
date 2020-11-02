import argparse
from functools import partial
from os.path import abspath, isfile

import itertools
from multiprocessing import Pool

import numpy as np
from numpy.lib.recfunctions import unstructured_to_structured

from qcore import srf, utils
from qcore.geo import ll2gp_multi, gp2ll

FLOAT_DTYPE = np.single


def ll_srf_dist_3d(srf_points, point):
    """
    Return 3d distance between the points taken from an srf and a given (lon, lat, dep) tuple point.
    :param srf_points: A structured array with the structure (lon, lat, dep)
    :param point: A tuple containing the (lon, lat, dep) of the point being queried
    :return: The distances of each srf point the given point in km
    """

    lon1 = np.radians(srf_points["lon"])
    lat1 = np.radians(srf_points["lat"])
    lon2 = np.radians(point[0])
    lat2 = np.radians(point[1])
    v_dist = srf_points["dep"] - point[2]

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    h_dist = 6378.139 * 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return np.sqrt(h_dist ** 2 + v_dist ** 2)


def calculate_damage(
    i_j_k,
    srf_points,
    begin_width_taper_km,
    end_width_taper_km,
    begin_depth_taper_km,
    end_depth_taper_km,
    min_damage_velocity,
    grid_spacing,
    global_nx,
    global_ny,
    domain_centre_lat,
    domain_centre_lon,
    domain_bearing,
):
    """
    Calculates the damage multiplier for a given grid point.
    If the grid point is between the maximum and minimum depth and/or distance values the damage multiplier is linearly scaled
    If it is between both, then the maximum damage multiplier is taken
    :param i_j_k: An x, y, z tuple giving the location of the grid point being tested
    :param srf_points: A structured array of srf points, with the structure (lon, lat, dep)
    :param begin_width_taper_km: The maximum distance at which the maximum damage occurs
    :param end_width_taper_km: The minimum distance at which no damage occurs
    :param begin_depth_taper_km: The maximum depth at which the maximum damage occurs
    :param end_depth_taper_km: The minimumm depth at which no damage occurs
    :param min_damage_velocity: The minimum damage multiplier possible
    :param grid_spacing: The size of the velocity model grid used (in km)
    :param global_nx: The number of grid points in the x direction
    :param global_ny: The number of grid points in the z direction
    :param domain_centre_lat: The longitude of the centre point
    :param domain_centre_lon: The latitude of the centre point
    :param domain_bearing: The rotation of the velocity model domain
    :return: If the grid point is to be damaged a (x, y, z, damage_multiplier) tuple giving the location and damage, otherwise None
    """
    i, j, k = i_j_k
    lon, lat = gp2ll(
        i,
        j,
        domain_centre_lat,
        domain_centre_lon,
        domain_bearing,
        global_nx,
        global_ny,
        grid_spacing,
    )

    dep = (k - 0.5) * grid_spacing
    dist = np.min(ll_srf_dist_3d(srf_points, (lon, lat, dep)))

    ratio = None

    # If distance_ratio is less than 0, we are in the maximum damage region
    # If distance_ratio is in the range (0, 1], we are in the taper region
    # If distance_ratio is greater than 1, we are beyond the damage region
    distance_ratio = (dist - begin_width_taper_km) / (
        end_width_taper_km - begin_width_taper_km
    )

    # If depth_ratio is less than 0, we are in the maximum damage depths
    # If depth_ratio is in the range (0, 1], we are in the taper depths
    # If depth_ratio is greater than 1, we are beyond the damage depths
    depth_ratio = (dep - begin_depth_taper_km) / (
        end_depth_taper_km - begin_depth_taper_km
    )

    if distance_ratio <= 0 and depth_ratio <= 0:
        ratio = min_damage_velocity
    elif distance_ratio <= 0 and depth_ratio <= 1:
        ratio = min_damage_velocity + depth_ratio * (1 - min_damage_velocity)
    elif distance_ratio <= 1 and depth_ratio <= 0:
        ratio = min_damage_velocity + distance_ratio * (1 - min_damage_velocity)
    elif distance_ratio <= 1 and depth_ratio <= 1:
        scale_ratio = max(
            depth_ratio,
            distance_ratio,
        )
        ratio = min_damage_velocity + scale_ratio * (1 - min_damage_velocity)
    if ratio is not None:
        return i, j, k, ratio


def check_points(
    bounding_box,
    srf_points,
    width_km,
    max_width_km,
    depth_km,
    max_depth_km,
    min_damage_velocity,
    vm_params,
    n=1,
):
    """
    For each VM point within the bounding box, check to see if any damage or velocity tapering needs to be applied
    :param bounding_box: a 3x2 tuple containing the ((min_nx, max_nx), (min_ny, max_ny), (min_nz, max_nz)) gridpoint bounds
    :param srf_points: The points in the srf as a structured array with the fields (lon, lat, dep)
    :param width_km: The maximum distance at which full damage is applied
    :param max_width_km: The minimum distance at which no damage is applied. Values in between are scaled linearly
    :param depth_km: The maximum depth at which full damage is applied
    :param max_depth_km: The mimimum depth at which no damage is applied. Values in between are scaled linearly
    :param min_damage_velocity: The minimum velocity multiplier
    :param vm_params: The dictionary loaded from the vm_params file
    :param n: The number of processes to use. Defaults to 1
    :return: A list of (x, y, z, ratio) tuples with damage to be applied
    """

    ((min_nx, max_nx), (min_ny, max_ny), (min_nz, max_nz)) = bounding_box

    global_nx = vm_params["nx"]
    global_ny = vm_params["ny"]
    grid_spacing = vm_params["hh"]
    domain_bearing = vm_params["MODEL_ROT"]
    domain_centre_lon = vm_params["MODEL_LON"]
    domain_centre_lat = vm_params["MODEL_LAT"]

    # Most of the input parameters are the same, the only difference being the VM point being queried
    func = partial(
        calculate_damage,
        srf_points=srf_points,
        begin_width_taper_km=width_km,
        end_width_taper_km=max_width_km,
        begin_depth_taper_km=depth_km,
        end_depth_taper_km=max_depth_km,
        min_damage_velocity=min_damage_velocity,
        grid_spacing=grid_spacing,
        global_nx=global_nx,
        global_ny=global_ny,
        domain_centre_lat=domain_centre_lat,
        domain_centre_lon=domain_centre_lon,
        domain_bearing=domain_bearing,
    )

    results = Pool(n).map(
        func,
        itertools.product(
            range(min_nx, max_nx), range(min_ny, max_ny), range(min_nz, max_nz)
        ),
    )
    return [x for x in results if x is not None]


def get_bounding_box(srf_corners, vm_params, max_width_km, max_depth_km):
    """
    Gets the outer bounds of the srf with respect to the velocity model grid, allowing most points too far from the srf to be filtered out
    :param srf_corners: The srf corners from srf.get_bounds
    :param vm_params: The dictionary loaded from the vm_params.yaml file
    :param max_width_km: The maximum damage width from the srf
    :param max_depth_km: The maximum damage depth from the surface
    :return: A 3x2 tuple containing the min, max - x, y, z bounds respectively
    """

    global_nx = vm_params["nx"]
    global_ny = vm_params["ny"]
    grid_spacing = vm_params["hh"]
    domain_bearing = vm_params["MODEL_ROT"]
    domain_centre_lon = vm_params["MODEL_LON"]
    domain_centre_lat = vm_params["MODEL_LAT"]

    # Get bounding box to filter out most points
    bounds_as_xy = ll2gp_multi(
        srf_corners,
        domain_centre_lon,
        domain_centre_lat,
        domain_bearing,
        global_nx,
        global_ny,
        grid_spacing,
    )

    min_nx, min_ny = np.min(np.floor(bounds_as_xy), axis=0)
    max_nx, max_ny = np.max(np.ceil(bounds_as_xy), axis=0)

    buffer = int(np.ceil(max_width_km / grid_spacing))

    min_nx = max(int(min_nx - buffer), 0)
    max_nx = min(int(max_nx + buffer), global_nx)
    min_ny = max(int(min_ny - buffer), 0)
    max_ny = min(int(max_ny + buffer), global_ny)
    min_nz, max_nz = 0, int(np.ceil(max_depth_km / grid_spacing))

    return (min_nx, max_nx), (min_ny, max_ny), (min_nz, max_nz)


def modify_file(results, pert_f_location, vm_params):
    """
    Writes the damage to the file, one result at a time

    :param results: a list of (x, y, z, damage) tuples to be applied to the file
    :param pert_f_location: The location of the perturbation file (must exist)
    :param vm_params: The dictionary loaded from the vm_params.yaml file
    """
    global_nx = vm_params["nx"]
    global_nz = vm_params["nz"]
    float_size = np.dtype(FLOAT_DTYPE).itemsize

    def locate_grid_point(x, y, z):
        """Converts x, y, z cordinates into the address"""
        return (global_nx * global_nz * y + global_nx * z + x) * float_size

    with open(pert_f_location, mode="br+") as pert_file:
        for i, j, k, val in results:
            location = locate_grid_point(i, j, k)
            pert_file.seek(location)
            existing_value = np.fromfile(pert_file, dtype=FLOAT_DTYPE, count=1)
            new_value: np.ndarray = val * existing_value
            pert_file.seek(-float_size, 1)
            new_value.tofile(pert_file)


def apply_fault_damage_zone(
    srf_location,
    vm_params,
    pert_f_location,
    depth_km,
    max_depth_km,
    width_km,
    max_width_km,
    min_damage_velocity,
    n_processes,
):
    """
    Applies the fault damage zone to the given perturbation file
    :param srf_location: The path of the srf
    :param vm_params: The dictionary loaded from the vm_params.yaml file
    :param pert_f_location: The location of the perturbation file
    :param depth_km: The maximum depth at which the maximum damage occurs
    :param max_depth_km: The minimum depth at which no damage occurs
    :param width_km: The maximum distance from the srf at which the maximum damage occurs
    :param max_width_km: The minimum distance from the srf at which no damage occurs
    :param min_damage_velocity: The minimum velocity multiplier
    :param n_processes: The number or processes to use for parallel code
    """
    # Get the bounding box to filter down the number of locations to check
    srf_corners = list(itertools.chain.from_iterable(srf.get_bounds(srf_location)))
    bounding_box = get_bounding_box(srf_corners, vm_params, max_width_km, max_depth_km)
    # Work out which ones to modify
    srf_points = unstructured_to_structured(
        srf.read_srf_points(srf_location),
        dtype=(np.dtype([("lon", float), ("lat", float), ("dep", float)])),
    )
    points_to_modify = check_points(
        bounding_box,
        srf_points=srf_points,
        width_km=width_km,
        max_width_km=max_width_km,
        depth_km=depth_km,
        max_depth_km=max_depth_km,
        min_damage_velocity=min_damage_velocity,
        vm_params=vm_params,
        n=n_processes,
    )
    # Apply the modification
    modify_file(points_to_modify, pert_f_location, vm_params)


def load_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("srf_location", type=abspath)
    parser.add_argument("vm_params_location", type=abspath)
    parser.add_argument("perturbation_file", type=abspath)
    parser.add_argument("-n", "--processes", default=1, type=int)

    add_fault_damage_zone_properties(parser)

    args = parser.parse_args()

    return args


def add_fault_damage_zone_properties(parser):
    parser.add_argument("--width_km", default=0.225, help="The core damage area", type=float)
    parser.add_argument(
        "--max_width_km",
        default=0.75,
        help="The distance to the outer edge of the taper region where no damage occurs",
        type=float,
    )
    parser.add_argument(
        "--max_velocity_drop",
        default=0.7,
        help="The maximum velocity multiplier. ie default = 0.7 results in 70% of Vs within the max damage zone",
        type=float,
    )
    parser.add_argument(
        "--depth_km",
        default=4.0,
        help="The maximum depth at which to apply full damage",
        type=float,
    )
    parser.add_argument(
        "--max_depth_km", default=6.0, help="The depth at which no damage occurs"
    )


def main():
    args = load_args()

    pert_f_location = args.perturbation_file

    vm_params = utils.load_yaml(args.vm_params_location)

    if not isfile(pert_f_location):
        # If the perturbation file doesn't exist yet, then we should make a new one of all ones (no perturbation)
        create_empty_perturbation_file(pert_f_location, vm_params)

    apply_fault_damage_zone(
        args.srf_location,
        vm_params,
        pert_f_location,
        args.depth_km,
        args.max_depth_km,
        args.width_km,
        args.max_width_km,
        args.max_velocity_drop,
        args.processes,
    )


def create_empty_perturbation_file(pert_f_location, vm_params):
    np.ones(
        shape=vm_params["nx"] * vm_params["ny"] * vm_params["nz"], dtype=FLOAT_DTYPE
    ).tofile(pert_f_location)


if __name__ == "__main__":
    main()

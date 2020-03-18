"""
Script used to generate points on a non-uniform grid of NZ based on population and vs30
"""

import argparse
import datetime
import os

import numpy as np
import pandas as pd

from qcore import formats, geo

BASE_GRID_SPACING = 8

VERSION = 2  # Increment this when creating a new production version
"""
Version 0 - 17.6
Version 1 - 18.6
Version 2 - 20.3
"""

pop_weight = 0.5
vs30_weight = 0.5
# At each grid size what does the score have to exceed to add new points to the grid
thresholds = [0.5, 0.67, 0.74, 0.77, 0.882, 0.9]
#             4km,  2km, 1km,  0.5,  0.25, 0.125


def score_vs30_delta(
    vs30_delta_value, vs30_delta_min_acceptable=0, vs30_delta_max_acceptable=600
):
    """
    :param vs30_delta_value: vs30 value to be scored
    :param vs30_delta_min_acceptable: minimum vs30 value in the scoring range (saturates at 0)
    :param vs30_delta_max_acceptable: maximum vs30 value in the scoring range (saturates at 1)
    :return: a score value between 0 and 1
    """
    if vs30_delta_value < vs30_delta_min_acceptable:
        score_vs30_delta = 0
    elif vs30_delta_value > vs30_delta_max_acceptable:
        score_vs30_delta = 1
    else:
        score_vs30_delta = (vs30_delta_value - vs30_delta_min_acceptable) / (
            vs30_delta_max_acceptable - vs30_delta_min_acceptable
        )
    return score_vs30_delta


def score_pop(
    population_value, min_population_acceptable=10, max_population_acceptable=30
):
    """
    :param population_value: population value to be scored
    :param min_population_acceptable: minimum population value in the scoring range (saturates at 0)
    :param max_population_acceptable: maximum population value in the scoring range (saturates at 1)
    :return: a score value between 0 and 1
    """
    if population_value > max_population_acceptable:
        score_population = 1.0
    elif population_value < min_population_acceptable:
        score_population = 0.0
    else:  # linear scaling between 0 and 1 for between min and max acceptable
        score_population = (population_value - min_population_acceptable) / (
            max_population_acceptable - min_population_acceptable
        )
    return score_population


def pop_distance(lat, lon, pop_df):
    """
    :param lat:
    :param lon:
    :param pop_df:
    :return: population / distance value for whole NZ.
    """
    locations = np.array([pop_df.lon.values, pop_df.lat.values]).T
    distances = geo.get_distances(locations, lon, lat)

    return (pop_df["pop"] / distances).sum()


def write_interim_data(outpath, lat, lon, **kwargs):
    """
    Logs scoring / value data to a folder for interogation / plotting
    :param outpath: Folder to contain files with data
    :param lat: latitude
    :param lon: longitude
    :param kwargs: Keyword argument pairs to log the values - each keyword is saved to a separate file
    """
    for key in kwargs:
        with open(os.path.join(outpath, key + ".csv"), "a") as fp:
            if not np.isnan(kwargs[key]):
                fp.write(f"{lon} {lat} {kwargs[key]}\n")


def write_ll_file(out_fname, points):
    """
    Writes out the points determined
    :param out_fname: Filename to write the non-uniform grid to
    :param points: dictionary containing the levels to be written
    """
    with open(out_fname, "w") as fp:
        for key in points:
            lvl_str = str(key) + str(VERSION)
            for i, (lat, lon) in enumerate(points[key]):
                name = lvl_str + (f"{i:x}").zfill(7 - len(lvl_str))
                fp.write(f"{lon} {lat} {name}\n")


def get_prospective_points(distance, lat, lon):
    """
    :param distance: The distance to create points in a NESW configuration from the lat/lon specified
    :param lat:
    :param lon:
    :return: a list of 4 points
    """
    prospective_points = [[] for i in range(4)]
    prospective_points[0] = geo.ll_shift(lat, lon, distance, 0)
    prospective_points[1] = geo.ll_shift(lat, lon, distance, 90)
    prospective_points[2] = geo.ll_shift(lat, lon, distance, 180)
    prospective_points[3] = geo.ll_shift(lat, lon, distance, 270)
    return prospective_points


def compute_point_score(lat, lon, prospective_points, pop_df, vs30_df, input_plots_dir=None):
    """
    Computes a score based on the population and vs30 datasets
    :param lat: point considered to be split
    :param lon: point considered to be split
    :param prospective_points: list of 4 points
    :param pop_df: dataframe containing population
    :param vs30_df: dataframe containing vs30
    :param input_plots_dir: folder to write input calculation data
    :return: score from 0 to 1
    """
    max_lat = prospective_points[0][0]
    min_lat = prospective_points[2][0]
    max_lon = prospective_points[1][1]
    min_lon = prospective_points[3][1]

    vs30_range = vs30_df[vs30_df.lon < max_lon][vs30_df.lon > min_lon][
        vs30_df.lat < max_lat
    ][vs30_df.lat > min_lat]
    vs30_delta = vs30_range.vs30.max() - vs30_range.vs30.min()
    vs30_delta_score = score_vs30_delta(vs30_delta)

    pop_dist = pop_distance(lat, lon, pop_df)
    pop_dist_score = score_pop(
        pop_dist, min_population_acceptable=2000, max_population_acceptable=9000
    )
    overall_score = pop_dist_score * pop_weight + vs30_delta_score * vs30_weight

    if input_plots_dir is not None:
        write_interim_data(
            input_plots_dir,
            lat,
            lon,
            vs30_delta=vs30_delta,
            pop_dist=pop_dist,
            vs30_delta_score=vs30_delta_score,
            pop_dist_score=pop_dist_score,
            overall_score=overall_score,
        )
    return overall_score


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("vs30_file", help="CSV containing lon/lat/vs30, no header")
    parser.add_argument("pop_file", help="CSV containing lon/lat/population, no header")
    parser.add_argument(
        "basegrid", help="SSV containing lon/lat/station_name, no header"
    )
    parser.add_argument(
        "output_filename", help="File to write to the lon/lat/station_name points to"
    )
    parser.add_argument("--hh", help="Minimum spacing between points", type=float, default="0.4")
    parser.add_argument("--input_plot_dir", help="directory to write output plots to")
    args = parser.parse_args()
    return args


def main(args):

    out_fname = args.output_filename
    hh = args.hh
    input_plots_dir = args.input_plot_dir
    if input_plots_dir is not None:
        input_plots_dir += str(datetime.datetime.now()).replace(" ", "_")
        os.makedirs(input_plots_dir)

    level = 0
    distance = BASE_GRID_SPACING

    vs30_df = pd.read_csv(
        args.vs30_file, header=None, names=["lon", "lat", "vs30"]
    ).dropna()
    base_grid = pd.read_csv(
        args.basegrid, header=None, names=["lon", "lat", "name"], delim_whitespace=True
    )
    pop_df = pd.read_csv(args.pop_file, names=["lon", "lat", "pop"])

    total_points = len(base_grid)
    print(f"{total_points} points in the basegrid")

    points = {}
    points[0] = list(zip(base_grid.lat, base_grid.lon))

    while distance / 2 > hh:
        distance /= 2
        prior_level = level
        threshold = thresholds[level]
        level += 1
        points[level] = []

        for lat, lon in points[prior_level]:
            prospective_points = get_prospective_points(distance, lat, lon)

            overall_score = compute_point_score(lat, lon, prospective_points, pop_df, vs30_df, input_plots_dir)

            if overall_score > threshold:
                for pt in prospective_points:
                    points[level].append(pt)

        total_points += len(points[level])
        print(f"{len(points[level])} points added to level {level}")

    write_ll_file(out_fname, points)
    print(f"\n{total_points} on non-uniform grid.")


if __name__ == "__main__":
    args = parse_args()
    main(args)

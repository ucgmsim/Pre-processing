import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('stat_file', help="output_statfile")
parser.add_argument('score_file', help="input score")

args = parser.parse_args()

thresholds = [0.4, 0.55, 0.65, 0.7, 0.882, 0.9]

score_df = pd.read_csv(args.score_file, delim_whitespace=True, names=['lon', 'lat', 'value'])
names_df = pd.read_csv('non_uniform_whole_nz-hh400.ll', delim_whitespace=True, names=['lon', 'lat', 'name'])

evaluation_df = score_df.merge(names_df)

points_in_basegrid = len(evaluation_df[[stat[0] == "0" for stat in evaluation_df.name.values]])
print(f"points in basegrid: {points_in_basegrid}")
total_points = points_in_basegrid

for i, threshold in enumerate(thresholds):
    stations_in_level = evaluation_df[[stat[0] == str(i) for stat in evaluation_df.name.values]]
    estimated_points = (stations_in_level.value >= threshold).sum() * 4
    total_points += estimated_points
    print(f"{len(stations_in_level)} points at level {i}")
    print(f"estimated {estimated_points} points meeting a {threshold} threshold at level {i}")

print(f"\ntotal estimated points: {total_points}")
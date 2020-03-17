import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("grid", help="Non uniform grid file", nargs='+')
args = parser.parse_args()

for grid in args.grid:
    stats = pd.read_csv(grid, delim_whitespace=True, header=None)

    results = {}
    total = 0

    for i, stat in stats.iterrows():
        level = stat[2][0]
        if level not in results:
            results[level] = 0
        results[level] += 1

    print(f"grid: {grid}")
    for key in results:
        print(f"Spacing: {8 / (2 ** int(key))}km points: {results[key]}")
        total += results[key]
    print(f"Total points: {total}")

#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("finest_level")
args = parser.parse_args()

input_file = args.input_file
finest_level = args.finest_level
with open(input_file) as f:
    lines = f.readlines()
    for line in lines:
        lon, lat, code = line.split()
        if int(code[0]) >= int(finest_level):
            continue
        print(f"{lon}\t{lat}\t{code.strip()}")

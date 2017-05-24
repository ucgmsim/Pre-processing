#!/usr/bin/env python

import sys


input_file = sys.argv[1]
finest_level = sys.argv[2]
with open(input_file) as f:
    lines = f.readlines()
    for line in lines:
        lon, lat, code = line.split("\t")
        if int(code[0]) >= int(finest_level):
            continue
        print "%s\t%s\t%s" %(lon, lat, code.strip())
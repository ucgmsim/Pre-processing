#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys
import tools
import generate_mesh
from tools import read_csv, mixed_criteria
import point_sorting
import hashlib
import ConfigParser
from StringIO import StringIO

usage = "%s [params_file] [centroids_file.csv] [vs500.csv]" % (sys.argv[0])
# Read all necessary csv files
if len(sys.argv) < 3:
    print sys.argv
    print "Please provide all the necessary input files"
    print usage
    exit(1)

# extract information from the parameters file
config = ConfigParser.ConfigParser()
with open(sys.argv[1]) as stream:
    stream = StringIO("[Parameters]\n" + stream.read())
    config.readfp(stream)

lon = float(config.get("Parameters","ORIGIN_LON").replace("'",""))
lat = float(config.get("Parameters","ORIGIN_LAT").replace("'",""))
rotate = float(config.get("Parameters","ORIGIN_ROT").replace("'",""))
x_len = float(config.get("Parameters","EXTENT_X").replace("'",""))
y_len = float(config.get("Parameters","EXTENT_Y").replace("'",""))
minimal_distance = float(config.get("Parameters","EXTENT_LATLON_SPACING").replace("'",""))
suffix = config.get("Parameters","SUFX").replace("'","")
maximal_distance = 8.0 # fixed to 8km

nx = int(x_len / maximal_distance) + 1
ny = int(y_len / maximal_distance) + 1
print nx, ny

dx = maximal_distance
dy = maximal_distance

cp = tools.ll2gp(lat, lon, lat, lon, rotate, nx, ny, maximal_distance)
print tools.gp2ll(cp[0], cp[1], lat, lon, rotate, nx, ny, maximal_distance)
print "Center point:", cp

initial_grid = generate_mesh.generate_regular_mesh(cp, nx, ny, dx, dy)
initial_domain = generate_mesh.compute_bounding_points(cp, nx, ny, dx, dy)
print "Uniform coarsest grid has %i points, %s" % (len(initial_grid), str(initial_domain))

parameters = {}
parameters["lon"] = lon
parameters["lat"] = lat
parameters["rot"] = rotate
parameters["nx"] = nx
parameters["ny"] = ny
parameters["dist"] = maximal_distance

population_information = {}
vs500_information = {}

population_information = tools.convert_csv_to_grid_points(sys.argv[2], parameters)
vs500_information = tools.convert_csv_to_grid_points(sys.argv[3], parameters, ordered=False)
vs_500_points = []
for coords, value in vs500_information.iteritems():
    vs_500_points.append([coords[0], coords[1], value])
vs500_sorted_grid = point_sorting.SortedPointGrid(initial_domain)
vs500_sorted_grid.subdivide(4)
# insert points
vs500_sorted_grid.insert_points(vs_500_points)

finest_regions = vs500_sorted_grid.return_finest_regions()
print len(finest_regions)
final_point_list = []
for region in finest_regions:
    points = region.points
    for p in points:
        final_point_list.append(p)

non_uniform_mesh = generate_mesh.NonUniformGrid(dx, minimal_distance, initial_grid, initial_domain)
print "Done Creating Grids"

weight_population = 0.5
weight_vs500 = 0.5
score_threshold = 0.5

f = lambda domain: mixed_criteria(domain, population_information, vs500_sorted_grid, weight_pop=weight_population,
                                  weight_vs=weight_vs500, threshold=score_threshold)

while non_uniform_mesh.refine(f):
    continue

print "Finished"
print non_uniform_mesh.current_distance

counter = 0
with open("non_uniform_%s%s_%s%s.ll" %(str(lon),str(lat),str(minimal_distance),suffix), "w") as fp:
    #fp.write("Lon,Lat,Level\n")
    for i in range(non_uniform_mesh.finest_level + 1):
        print "writing level", i
        levels = non_uniform_mesh.get_points_at_level(i)
        for p in levels:
            longitude, latitude = tools.gp2ll(p[0], p[1], lat, lon, rotate, nx, ny, maximal_distance)
            # fp.write("%f, %f, %f\n" % (p[0], p[1], float(i)))
            counter += 1
            # TODO: this only creates an 8 char name, make it better
            name = hashlib.md5("%f.5%f.5" % (longitude, latitude)).hexdigest()[:8]
            fp.write("%f\t%f\t%s\n" % (longitude, latitude, name ))


print "Non-uniform grid has %i points" % counter

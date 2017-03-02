import json
import sys
from shapely.geometry import Point, shape
import csv
import tools

# load GeoJSON file containing sectors
# with open('test.json',mode='r') as f:
#     js = json.load(f)

# load VeloModCorners.txt
list_lon = []
list_lat = []
with open("VeloModCorners.txt") as corners:
    lines = corners.readlines()
    for line in lines:
        try:
            lon, lat = map(float, line.split("\t"))
            print lon, lat
            list_lon.append(lon)
            list_lat.append(lat)
        except ValueError:
            print "Ignoring",line.strip()

geo_json = {"type": "Feature",
             "geometry": {
                 "type": "Polygon",
                 "coordinates" : [[[lon,lat] for lon,lat in zip(list_lon,list_lat)]]
             }
            }
print geo_json

point_list = []
all_points = tools.read_csv(sys.argv[1], delimiter="\t")
for point in all_points:
    point_list.append(Point(float(point[0]), float(point[1])))

print len(point_list)

with open('cropped_' + sys.argv[1], mode="w") as out:
    # check each polygon to see if it contains the point
    for point in point_list:
        polygon = shape(geo_json['geometry'])

        if polygon.contains(point):
           out.write("%f %f\n" % (point.x, point.y))

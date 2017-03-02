import json
import sys
from shapely.geometry import Point, shape
import csv
import tools

# load GeoJSON file containing sectors
with open('test.json',mode='r') as f:
    js = json.load(f)

point_list = []
all_points = tools.read_csv(sys.argv[1],delimiter="\t")
for point in all_points:
    point_list.append(Point(float(point[0]),float(point[1])))

print len(point_list)

with open('cropped_'+sys.argv[1],mode="w") as out:
    # check each polygon to see if it contains the point
    for point in point_list:
        for feature in js['features']:
            polygon = shape(feature['geometry'])
            if polygon.contains(point):
                out.write("%f %f\n" % (point.x, point.y))



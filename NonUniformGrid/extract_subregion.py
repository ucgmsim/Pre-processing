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
poly_file = sys.argv[2]
with open(poly_file) as corners:
    lines = corners.readlines()
    for line in lines:
        try:
            lat, lon = map(float, line.split(","))
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
#all_points = tools.read_csv(sys.argv[1], delimiter="\t")

with open(sys.argv[1]) as f:
    lines = f.readlines()
    lines = map(lambda x: x.strip().replace("\t"," "), lines)
    all_points = map(lambda x: x.split(" "), lines)

for point in all_points:
    point = filter(lambda x: x != "", point)
    print point
    point_list.append([float(point[0]), float(point[1]), point[2]])


print len(point_list)

with open('cropped_test.csv' , mode="w") as out:
    # check each polygon to see if it contains the point
    polygon = shape(geo_json['geometry'])
    for point in point_list:
        tmp_point = Point(point[0],point[1])
        if polygon.contains(tmp_point):
           out.write("%f %f %s\n" % (point[0], point[1], point[2]))

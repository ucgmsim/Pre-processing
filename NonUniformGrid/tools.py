import math
import csv

# model reference location
from subprocess import Popen, PIPE


def read_csv(filename, delimiter=","):
    with open(filename, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=delimiter)
        result = []
        for each_row in reader:
            result.append(each_row)
        return result


def dump_csv_from_point_list(csv_file_name, points, header=None):
    with open( csv_file_name, "w") as output_file:
        if header:
            output_file.write("%s\n" %header)

        for point in points:
            try:
                output_file.write(",".join(map(str,point))+"\n")
            except:
                print "error dumping",point


# Uses population csv and vs_500 csv
def mixed_criteria(domain, population_hash, vs500_sorted_grid, min_population_acceptable=10.0, weight_pop=0.5,
                   weight_vs=0.5, threshold=0.5):
    interest_points_population = {}
    interest_points_vs500 = {}
    # TODO: change to use sorted points
    for coords, value in population_hash.iteritems():
        if point_in_rectangle(coords[0], coords[1], domain):
            interest_points_population[coords] = value

    regions = vs500_sorted_grid.regions_around_domain(domain)
    for region in regions:
        points = region.points
        for point in points:
            if point_in_rectangle(point[0], point[1], domain):
                interest_points_vs500[coords] = point[2]

    # analyze all the points attributing score
    population_value = 0.0
    counter = 0
    # Population: mean population in the domain
    for point, value in interest_points_population.iteritems():
        if population_value < value:
            population_value = value

    vs500_value = 100.0
    # vs_500: min value in the domain
    for point, value in interest_points_vs500.iteritems():
        if value < vs500_value:
            vs500_value = value

    # vs_500: finest grid if 0.5 and coarsest if 2.0
    # population: finest grid if population > 2*min_population, coarsest if < min_population
    vs500_max_acceptable = 2.0
    vs500_min_acceptable = 0.5

    max_population_acceptable = 3.0 * min_population_acceptable
    if population_value > max_population_acceptable:
        score_population = 1.0
    elif population_value < min_population_acceptable:
        score_population = 0.0
    else:
        score_population = 1.0 - (
        (max_population_acceptable - population_value) / (max_population_acceptable - min_population_acceptable))

    if vs500_value > vs500_max_acceptable:
        score_vs500 = 0.0
    elif vs500_value < vs500_min_acceptable:
        score_vs500 = 1.0
    else:
        score_vs500 = (vs500_max_acceptable - vs500_value) / (vs500_max_acceptable - vs500_min_acceptable)

    total_score = weight_pop * score_population + weight_vs * score_vs500
    assert(total_score <= 1.0)

    if total_score >= threshold:
        return True
    return False



def convert_csv_to_grid_points(csv_file_name, parameters, ordered=True):
    """
    Convert a csv file of the form Lon,Lat,Property into grid points.
    :param csv_file_name: file name of the csv
    :param parameters: a hash containing all the parameters needed for the ll2gp conversion
    :param ordered: True if each line is Lon,Lat,Prop, otherwise if it is Lat,Lon,Prop should be false
    :return:
    """
    lines = read_csv(csv_file_name)
    lat_input = []
    lon_input = []
    property_input = []
    result = {}
    for line in lines:
        try:
            if ordered:
                lon_point, lat_point, prop = line
            else:
                lat_point, lon_point, prop = line

            property_input.append(float(prop))
            lat_input.append(lat_point)
            lon_input.append(lon_point)

        except:
            continue

    x_coords, y_coords = ll2gp_file(lat_input, lon_input, parameters["lat"], parameters["lon"], parameters["rot"], \
                                    parameters["nx"], parameters["ny"], parameters["dist"])

    for i in range(len(x_coords)):
        try:
            if x_coords[i] and y_coords[i]:
                result[x_coords[i], y_coords[i]] = property_input[i]
        except Exception:
            print Exception.message

    return result


def point_in_rectangle(px, py, rectangle):
    """
    Check if a point is included into a regular rectangle that aligned with the x-y axis
    :param px: x value of the point
    :param py: y value of the point
    :param rectangle: the coordinates of the 4 corners of the rectangle
    :return: True or False depending if the point is in poly or not
    """
    max_x_val = float("-inf")
    max_y_val = float("-inf")
    min_x_val = float("inf")
    min_y_val = float("inf")

    assert (len(rectangle) == 4)

    for p in rectangle:
        if p[0] >= max_x_val:
            max_x_val = p[0]
        if p[1] >= max_y_val:
            max_y_val = p[1]
        if p[0] <= min_x_val:
            min_x_val = p[0]
        if p[1] <= min_y_val:
            min_y_val = p[1]

    return px >= min_x_val and px <= max_x_val and py >= min_y_val and py <= max_y_val


def ll_dist(lon1, lat1, lon2, lat2):
    """
    Return distance between a pair of lat, lon points.
    """
    R_EARTH = 6378.139
    lat1, lat2, dlon, dlat = map(math.radians, \
                                 [lat1, lat2, (lon2 - lon1), (lat2 - lat1)])

    a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    return R_EARTH * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def ll_dist(p0, p1):
    return ll_dist(p0[0], p0[1], p1[0], p1[1])


def ll2gp_file(lat_vals, lon_vals, mlat, mlon, rot, nx, ny, hh):
    # derived parameters
    xlen = (nx - 1) * hh
    ylen = (ny - 1) * hh
    # where the x plane points
    xazim = (rot + 90) % 360

    assert (len(lat_vals) == len(lon_vals))

    in_file_name = "in.tmp"
    with open(in_file_name, "w") as fp:
        for i in range(len(lat_vals)):
            fp.write("%s %s\n" % (lon_vals[i], lat_vals[i]))

    # run binary, get output
    # output is displacement (x, y) from center (mlon, mlat), in kilometres
    cmd = ["./ll2xy", 'mlat=%s' % (mlat), 'mlon=%s' % (mlon), \
           'geoproj=1', 'center_origin=1', 'xazim=%s' % (xazim), 'xlen=%s' % (xlen), 'ylen=%s' % (ylen),
           "infile=%s" % (in_file_name)]
    print ' '.join(cmd)
    p_conv = Popen(cmd,
                   stdin=PIPE, stdout=PIPE)
    # stdout = p_conv.communicate('%s %s' % (lon, lat))[0]
    # x, y = map(float, stdout.split())

    points_x = []
    points_y = []

    max_x = (nx - 1) * hh
    max_y = (ny - 1) * hh

    for line in p_conv.stdout:
        try:
            x, y = filter(lambda x: x != "", line.split(" "))
            x, y = map(float, [x, y])
        except:
            continue

        x += (max_x) * 0.5
        y += (max_y) * 0.5
        # then convert back to grid spacing
        # TODO: this is dodgy for my case
        # x /= hh
        # y /= hh
        # gridpoints are discrete
        x = int(round(x))
        y = int(round(y))

        # nx values range from 0 -> nx - 1
        if not (-1 < x < xlen) or not (-1 < y < ylen):
            # raise InputError('Input outside simulation domain.')
            # print 'Warning: Input outside simulation domain', x, y
            points_x.append(None)
            points_y.append(None)
            continue

        points_x.append(x)
        points_y.append(y)

    return points_x, points_y


def ll2gp(lat, lon, mlat, mlon, rot, nx, ny, hh):
    """
    Converts latitude/longitude to a gridpoint position.
    Three main modes of operation:
    1: No dx, dy (= 1): gridpoint
    2: dx or dy != 1: closest gridpoint considering dx/dy
    3: decimated: gridpoint number if only decimated points existed
    """
    # derived parameters
    xlen = nx * hh
    ylen = ny * hh
    # where the x plane points
    xazim = (rot + 90) % 360

    # run binary, get output
    # output is displacement (x, y) from center, in kilometres
    cmd = ["./ll2xy", 'mlat=%s' % (mlat), 'mlon=%s' % (mlon), \
           'geoproj=1', 'center_origin=1', 'h=%s' % (hh), \
           'xazim=%s' % (xazim), 'xlen=%s' % (xlen), 'ylen=%s' % (ylen)]
    print ' '.join(cmd)
    p_conv = Popen(cmd,
                   stdin=PIPE, stdout=PIPE)
    stdout = p_conv.communicate('%s %s' % (lon, lat))[0]
    x, y = map(float, stdout.split())

    # convert displacement to grid points
    # first make the distance relative to top corner
    # has to be 'nx - 1', it's just how ll2xy works.
    # nx = 1400 means the largest value is 1399.
    max_x = (nx - 1) * hh
    max_y = (ny - 1) * hh
    x += (max_x) * 0.5
    y += (max_y) * 0.5

    # gridpoints are discrete
    x = int(round(x))
    y = int(round(y))
    print x, y
    # nx values range from 0 -> nx - 1
    if not (-1 < x < xlen) or not (-1 < y < ylen):
        # raise InputError('Input outside simulation domain.')
        print ('Warning: Input outside simulation domain')
        return None

    return x, y


def gp2ll(x, y, mlat, mlon, rot, nx, ny, hh):
    """
    Converts a gridpoint position to latitude/longitude.
    """
    # derived parameters
    xlen = nx * hh
    ylen = ny * hh
    # where the x plane points
    xazim = (rot + 90) % 360

    # convert gridpoint to offset km
    max_x = (nx - 1) * hh
    max_y = (ny - 1) * hh
    # length from corner
    # x *= hh
    # y *= hh
    # length from centre origin
    x -= max_x * 0.5
    y -= max_y * 0.5

    # run binary, get output
    cmd = ['./xy2ll', 'mlat=%s' % (mlat), 'mlon=%s' % (mlon), \
           'geoproj=1', 'center_origin=1', 'h=%s' % (hh), \
           'xazim=%s' % (xazim), 'xlen=%s' % (xlen), 'ylen=%s' % (ylen)]

    p_conv = Popen(cmd, stdin=PIPE, stdout=PIPE)
    stdout = p_conv.communicate('%s %s' % (x, y))[0]
    # lon, lat
    return map(float, stdout.split())

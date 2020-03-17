from tools import point_in_rectangle
from multiprocessing import pool


def test_point_to_refine(packed_values):
    point, next_distance, criteria = packed_values
    print point, next_distance, criteria
    tentative_points = [(point[0] - next_distance, point[1] - next_distance), \
                        (point[0] + next_distance, point[1] - next_distance), \
                        (point[0] + next_distance, point[1] + next_distance), \
                        (point[0] - next_distance, point[1] + next_distance),
                        (point[0] - next_distance, point[1]), \
                        (point[0], point[1] - next_distance), \
                        (point[0] + next_distance, point[1]), \
                        (point[0], point[1] + next_distance)]
    tentative_domain = [(point[0] - next_distance, point[1] - next_distance), \
                        (point[0] + next_distance, point[1] - next_distance), \
                        (point[0] + next_distance, point[1] + next_distance), \
                        (point[0] - next_distance, point[1] + next_distance)]

    domain_accepted = criteria(tentative_domain)

    if domain_accepted:
        for tentative_point in tentative_points:
            # ignore potential outside points
            if point_in_rectangle(tentative_point[0], tentative_point[1], self.domain):
                return tentative_point
    return []


def generate_regular_mesh(center_point, nx, ny, distance_x, distance_y):
    points = []
    nb_points_x = nx
    nb_points_y = ny

    if nb_points_x - int(nb_points_x) > 0.0:
        print nb_points_x, " not exact"
    if nb_points_y - int(nb_points_y) > 0.0:
        print nb_points_y, " not exact"

    for i in range(nb_points_x):
        for j in range(nb_points_y):
            points.append([center_point[0] + (i + 0.5) * distance_x - (nx * distance_x) / 2.0,
                           center_point[1] + (j + 0.5) * distance_y - (ny * distance_y) / 2.0])

    return points


def compute_bounding_points(center_point, nx, ny, distance_x, distance_y):
    nb_points_x = nx
    nb_points_y = ny
    if nb_points_x - int(nb_points_x) > 0.0:
        print nb_points_x, " not exact"
    if nb_points_y - int(nb_points_y) > 0.0:
        print nb_points_y, " not exact"

    return [[center_point[0] + (0.5) * distance_x - (nx * distance_x) / 2.0,
             center_point[1] + (0.5) * distance_y - (ny * distance_y) / 2.0],
            [center_point[0] + (nx + 0.5) * distance_x - (nx * distance_x) / 2.0,
             center_point[1] + (0.5) * distance_y - (ny * distance_y) / 2.0],
            [center_point[0] + (nx + 0.5) * distance_x - (nx * distance_x) / 2.0,
             center_point[1] + (ny + 0.5) * distance_y - (ny * distance_y) / 2.0],
            [center_point[0] + (0.5) * distance_x - (nx * distance_x) / 2.0,
             center_point[1] + (ny + 0.5) * distance_y - (ny * distance_y) / 2.0]
            ]


class NonUniformGrid:
    levels = {}

    def __init__(self, initial_distance, min_distance, initial_grid, initial_domain):
        self.current_distance = initial_distance
        self.min_distance = min_distance

        self.finest_level = 0
        self.levels[self.finest_level] = initial_grid
        self.domain = initial_domain

    def refine_multi(self, criteria):
        current_finest_level = self.finest_level
        print "refining level", current_finest_level
        next_distance = self.current_distance / 2.0

        if next_distance < self.min_distance:
            print next_distance, "vs.", self.min_distance
            return False
        p = pool.Pool(4)
        current_points = self.levels[current_finest_level]
        packed_values = []
        for point in current_points:
            packed_values.append([point, next_distance, criteria])
        new_mesh = p.map(test_point_to_refine, packed_values)
        new_mesh = list(set(new_mesh))
        print "Added: ", len(new_mesh), "new points to level", self.finest_level
        self.levels[self.finest_level] = new_mesh
        return True

    def refine(self, criteria, threshold):
        current_finest_level = self.finest_level
        print "refining level", current_finest_level
        next_distance = self.current_distance / 2.0

        if next_distance < self.min_distance:
            print next_distance, "vs.", self.min_distance
            return False
        new_mesh = []
        counter = 0
        for point in self.levels[current_finest_level]:
            # compute 4 new points
            counter += 1
            tentative_points = [(point[0] - next_distance, point[1] - next_distance), \
                                (point[0] + next_distance, point[1] - next_distance), \
                                (point[0] + next_distance, point[1] + next_distance), \
                                (point[0] - next_distance, point[1] + next_distance),
                                (point[0] - next_distance, point[1]), \
                                (point[0], point[1] - next_distance), \
                                (point[0] + next_distance, point[1]), \
                                (point[0], point[1] + next_distance)]
            tentative_domain = [(point[0] - next_distance, point[1] - next_distance), \
                                (point[0] + next_distance, point[1] - next_distance), \
                                (point[0] + next_distance, point[1] + next_distance), \
                                (point[0] - next_distance, point[1] + next_distance)]

            domain_accepted = criteria(tentative_domain, threshold)

            if domain_accepted:
                for tentative_point in tentative_points:
                    # ignore potential outside points
                    if point_in_rectangle(tentative_point[0], tentative_point[1], self.domain):
                        new_mesh.append(tentative_point)
            if counter % 100 == 0:
                print "Points dealt with", counter

        self.current_distance = next_distance
        self.finest_level += 1
        unique_points = set(new_mesh)
        new_mesh = list(unique_points)
        print "Added: ", len(new_mesh), "new points to level", self.finest_level
        self.levels[self.finest_level] = new_mesh
        return True

    def get_points_at_level(self, level):
        # type: (int) -> list
        """

        :type level: list
        """
        if level > self.finest_level:
            return None
        return self.levels[level]

    def get_all_points(self):
        result = []
        for level in self.levels.items():
            for point in level:
                result.append(point)

        return result

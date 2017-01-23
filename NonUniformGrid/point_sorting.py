
from tools import point_in_rectangle


def compute_middle_point(p0, p1):
    """
    Compute the middle point between 2 points by simple arithmetical mean
    :param p0: first point as a list or tuple
    :param p1: second point
    :return: a list that contains the middle point
    """
    return [(p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0]


def subdivide_domain(domain):
    assert (len(domain) == 4)
    p0, p1, p2, p3 = domain
    p01 = compute_middle_point(p0, p1)
    p03 = compute_middle_point(p0, p3)
    p12 = compute_middle_point(p1, p2)
    p23 = compute_middle_point(p2, p3)
    # special case, find gravity center of 4 points
    p02 = compute_middle_point(p03, p12)

    domain1 = [p0, p01, p02, p03]
    domain2 = [p01, p1, p12, p02]
    domain3 = [p02, p12, p2, p23]
    domain4 = [p03, p02, p23, p3]
    return [domain1, domain2, domain3, domain4]


def domains_intersect(domain1, domain2):
    for point in domain1:
        if point_in_rectangle(point[0],point[1],domain2):
            return True
    for point in domain2:
        if point_in_rectangle(point[0],point[1],domain1):
            return True

    return False


class GridRegion:
    def __init__(self, domain):
         self.domain = domain
         self.points = []

    def insert_point(self, point):
        self.points.append(point)

    def __str__(self):
        return "domain: "+str(self.domain[0])+","+str(self.domain[1])+","+str(self.domain[2])+","+str(self.domain[3])

    @property
    def domain(self):
        return self.domain

    @property
    def points(self):
        return self.points


class SortedPointGrid:
    """

    """
    levels = []

    def __init__(self, original_domain):
        self.original_domain = original_domain
        self.levels.append([GridRegion(original_domain)])
        self.max_level = 0

    def subdivide(self, max_level):

        for current_level in range(1,max_level+1):
            self.levels.append([])
            for grid_region in self.levels[current_level-1]:
                # subdivide
                new_domains = subdivide_domain(grid_region.domain)
                # put new domains
                for new_domain in new_domains:
                    self.levels[current_level].append(GridRegion(new_domain))

        self.max_level = max_level

    def close_points(self, point):
        regions = self.levels[self.max_level]
        print "going through %i regions"%(len(regions))
        for r in regions:
            if point_in_rectangle(point[0],point[1],r.domain):
                print r.points,r.domain
                return r.points

        return None

    def regions_around_domain(self, domain):
        regions = self.levels[self.max_level]
        result = []
        #print "going through %i regions"%(len(regions))
        for r in regions:
            if domains_intersect(domain,r.domain):
                result.append(r)

        #print "Found %i regions" %(len(result))
        return result

    def return_finest_regions(self):
        return self.levels[self.max_level]

    def insert_points(self, points):
        # we always insert at max level
        total_points = len(points)
        inserted_points = 0
        missing_points = 0
        for point in points:
            point_inserted = False
            for grid_region in self.levels[self.max_level]:
                if point_in_rectangle(point[0],point[1],grid_region.domain):
                    grid_region.insert_point(point)
                    #print str(point),"inserted on ",grid_region.domain
                    inserted_points += 1
                    point_inserted = True
                    break

            if not point_inserted:
                print point,"not inserted"
                missing_points += 1
        print "Finished inserting points"
        assert(inserted_points + missing_points == total_points)

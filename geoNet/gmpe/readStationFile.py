# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 13:29:34 2016

@author: rmc84

A collection of functions relating to reading the station and rupture files, taken from the matlab code.
In each function, the name of the original .m function is indicated.
Function and variable names from matlab are mostly preserved.

"""

from math import radians, sin, cos, asin, sqrt, ceil
import numpy as np

class Points:
    def __init__(self):
         self.Lat = []
         self.Lon = []
         self.Depth = []
         self.Name = []

class Points_np(object):
    def __init__(self, Points):
        self.Lat = np.asarray(Points.Lat)*np.pi/180.
        self.Lon = np.asarray(Points.Lon)*np.pi/180.
        self.Depth = np.asarray(Points.Depth)
        self.Name = Points.Name

def horizdist(loc1, loc2_lat, loc2_lon):
    """From ComputeSourceToSiteDistance.m """
    # computes great circle distance between 2 set of (lat, lng) (in degrees)
    EARTH_RADIUS_MEAN = 6371.0072   # Authalic mean radius in km

    # calculations are all in radians
    lat_1, lon_1, lat_2, lon_2  = map(radians, (loc1.Lat[0], loc1.Lon[0], loc2_lat, loc2_lon))
    lat = lat_2 - lat_1
    lon = lon_2 - lon_1
    d = sin(lat * 0.5) ** 2 + cos(lat_1) * cos(lat_2) * sin(lon * 0.5) ** 2
    h = 2.0 * EARTH_RADIUS_MEAN * asin(sqrt(d))
    return h

def horizdist_np(loc1, loc2_lat, loc2_lon):
    
    EARTH_RADIUS_MEAN = 6371.0072   # Authalic mean radius in km

    # calculations are all in radians
    #everything is in radians
    #lat_1, lon_1, lat_2, lon_2  = map(radians, (loc1.Lat[0], loc1.Lon[0], loc2_lat, loc2_lon))

    lat_1, lon_1, lat_2, lon_2  = (loc1.Lat[0], loc1.Lon[0], loc2_lat, loc2_lon)
    lat = lat_2 - lat_1
    lon = lon_2 - lon_1
    d = np.sin(lat * 0.5) ** 2 + np.cos(lat_1) * np.cos(lat_2) * np.sin(lon * 0.5) ** 2
    h = 2.0 * EARTH_RADIUS_MEAN * np.arcsin(np.sqrt(d))
    return h

def readStationCordsFile(station_file):
    """Based on readStationCordsFile.m
    """
    try:
        fp = open(station_file, 'r')
    except IOError:
        print 'Station filename is not valid. Returning from function readStationCordsFile'
        raise
        return

    stations = Points()
    # python optimisation
    add_lat = (stations.Lat).append
    add_lon = (stations.Lon).append
    add_name = (stations.Name).append

    for line in fp:
        station_info = line.split()
        add_lon(float(station_info[0]))
        add_lat(float(station_info[1]))
        add_name(station_info[2])
    stations.Depth = [0.0] * len(stations.Lon)

    fp.close()
    return stations

def readSrfFile(srfFile):
    """Based on readSrfFile.m
    Read the SRF (Standard Rupture Format)
    The rupture velocities are not read in presently as not used
    """

    try:
        fp = open(srfFile, 'r')
    except IOError:
        print('SRF filename is not valid. Returning from function readSrfFile')
        raise
        return
    # reading line by line is slightly slower and next() doesn't work with readline()
    # won't use much memory for really large files however. ~0.1 seconds / 10.6 slower
    next_line = fp.readline

    finite_fault = Points()
    # python optimisation
    add_lat = (finite_fault.Lat).append
    add_lon = (finite_fault.Lon).append
    add_depth = (finite_fault.Depth).append

    # skip version line
    next_line()
    # number of double line entries in 'PLANE' section
    n_seg = int(next_line().split()[1])

    # skip 'PLANE' section, onto 'POINTS' section
    for _ in xrange(n_seg):
        next_line()
        next_line()

    # now loop through rest of the file and get lon, lat and depth
    for _ in xrange(int(next_line().split()[1])):
        # index 0-7: lon, lat, depth, stk, dip, area, tInit, dt
        info_line_1 = next_line().split()
        # index 0-6: rake, slip1, nt1, slip2, nt2, slip3, nt3
        info_line_2 = next_line().split()

        add_lon(float(info_line_1[0]))
        add_lat(float(info_line_1[1]))
        add_depth(float(info_line_1[2]))

        # skip values for point
        for _ in xrange(int(ceil(int(info_line_2[2])/6.0))):
            next_line()

    fp.close()
    return finite_fault

def computeSourcetoSiteDistance(FiniteFault, Site):
    """ Purpose: compute the distance in km from the finite fault plane to the
    site (of an instrument or other).
    Based on ComputeSourceToSiteDistance.m """

    # start values, no distance should be longer than this
    Rjb= 99999
    Rrup=99999

    # initialize the ith Fault
    faulti_lat = 0.0
    faulti_lon = 0.0
    faulti_depth = 0.0

    # for subfaults, calculate distance, update if shortest
    for i in range(len(FiniteFault.Lat)):
        faulti_lat = FiniteFault.Lat[i]
        faulti_lon = FiniteFault.Lon[i]
        faulti_depth = FiniteFault.Depth[i]

        h = horizdist(Site, faulti_lat, faulti_lon)
        v = Site.Depth[0] - faulti_depth

        if abs(h) < Rjb:
            Rjb = h

        d = sqrt(h ** 2 + v ** 2)
        if d < Rrup:
            Rrup = d

    return Rrup, Rjb


def computeSourcetoSiteDistance_np(FiniteFault, Site):
    """ Purpose: compute the distance in km from the finite fault plane to the
    site (of an instrument or other).
    Based on ComputeSourceToSiteDistance.m """

    # start values, no distance should be longer than this
    Rjb= 99999
    Rrup=99999

    # for subfaults, calculate distance, update if shortest
    h = horizdist_np(Site, FiniteFault.Lat, FiniteFault.Lon)
    v = Site.Depth[0] - FiniteFault.Depth

    #if abs(h) < Rjb:
    #    Rjb = h

    #assert np.all( h < Rjb)
    Rjb = h.min()
    d = np.sqrt(h*h + v*v)
    #if d < Rrup:
    #    Rrup = d
    Rrup = d.min()

    return Rrup, Rjb

def computeRrup(stationFile, srfFile):
    """Wrapper function to calculate the rupture distance from the station and srf files"""

    # read in list of stations
    try:
        stations = readStationCordsFile(stationFile)      #dodgy var name from matlab
    except IOError:
        print 'Station filename is not valid. Returning from function computeRrup'
        raise
        return

    # read in the rupture file
    try:
        FiniteFault = readSrfFile(srfFile)
    except IOError:
        print 'SRF filename is not valid. Returning from function computeRrup'
        raise
        return

    # initialize the station and fault class
    Site = Points()

    # Number of stations from the station file
    n_stations = len(stations.Lat)
    # initialize the rupture distances list
    r_rups = []
    r_jbs = []

    # loop over the stations
    for j in xrange(n_stations):
        try:
            Site.Lat[0] = stations.Lat[j]
            Site.Lon[0] = stations.Lon[j]
            Site.Depth[0] = stations.Depth[j]
        except IndexError:
            Site.Lat.append(stations.Lat[j])
            Site.Lon.append(stations.Lon[j])
            Site.Depth.append(stations.Depth[j])
        rrup_rjbs = computeSourcetoSiteDistance(FiniteFault,Site)
        r_rups.append(rrup_rjbs[0])
        r_jbs.append(rrup_rjbs[1])

    # return the source to site rupture distances
    return (r_rups, r_jbs, n_stations, stations.Name)



def computeRrup_np(stationFile, srfFile):
    """Wrapper function to calculate the rupture distance from the station and srf files"""

    # read in list of stations
    try:
        stations = readStationCordsFile(stationFile)      #dodgy var name from matlab
    except IOError:
        print 'Station filename is not valid. Returning from function computeRrup'
        raise
        return

    # read in the rupture file
    try:
        FiniteFault = readSrfFile(srfFile)
    except IOError:
        print 'SRF filename is not valid. Returning from function computeRrup'
        raise
        return

    #all changes from degrees to radians
    stations = Points_np(stations)
    FiniteFault = Points_np(FiniteFault)
    # initialize the station and fault class
    Site = Points()

    # Number of stations from the station file
    n_stations = len(stations.Lat)
    # initialize the rupture distances list
    r_rups = []
    r_jbs = []

    # loop over the stations
    for j in xrange(n_stations):
        try:
            Site.Lat[0] = stations.Lat[j]
            Site.Lon[0] = stations.Lon[j]
            Site.Depth[0] = stations.Depth[j]
        except IndexError:
            Site.Lat.append(stations.Lat[j])
            Site.Lon.append(stations.Lon[j])
            Site.Depth.append(stations.Depth[j])
        rrup_rjbs = computeSourcetoSiteDistance_np(FiniteFault,Site)
        r_rups.append(rrup_rjbs[0])
        r_jbs.append(rrup_rjbs[1])

    # return the source to site rupture distances
    return (r_rups, r_jbs, n_stations, stations.Name)

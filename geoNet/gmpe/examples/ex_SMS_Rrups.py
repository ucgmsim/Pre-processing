"""
Read stations_list.ll and find the Rrups for each SMS in it.
"""
import numpy as np
from geoNet.utils import read_statsll
from geoNet.gmpe import readStationFile as rsf
from geoNet.gmpe.calculateGMPE import set_faultprop, set_siteprop


# Read the standard rupture file
srf_fname = "/nesi/projects/nesi00213/RupModel/ahsan/2011Apr29_m4pt9/ptSource/Srf/2011Apr29_v2_m4pt9.srf"
FiniteFault = rsf.readSrfFile(srf_fname)
FiniteFault = rsf.Points_np(FiniteFault)


stats_floc = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017"
stats_fname = "20110429_190804_eventStats_2017-01-10.ll"
event_stats = read_statsll(stats_floc, stats_fname)

print("Station code: Rrup, Rjb")
for stat_code in event_stats:
    lon, lat = event_stats[stat_code]

    site = rsf.Points()
    site.Lon.append(lon)
    site.Lat.append(lat)
    site.Depth.append(0.0)
    site = rsf.Points_np(site)

    Rrup, Rjb = rsf.computeSourcetoSiteDistance_np(FiniteFault, site)

    print("{:^10s} {:^15.4f} {:^15.4f}".format(stat_code, Rrup, Rjb))

import os
from time import time
from geoNet.geoNet import scrapeGeoNet as sg
from geoNet.geoNet.gen_stats_kml import write_stats_kml

init_time=time()


#Firs get statsll dictionary by assigning lon, lat from the saved list known 
#to be in WGS84 coordinates
fname="ex_stats.ll"
loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo"
fname_all_geoNet_stats="all_geoNet_stats+2016-12-20.ll"
loc_V1A="/".join([os.getcwd(), "Vol1", "data"])
(
event_stats1,
fname_statsll
) = sg.statsll_from_V1A_v2(
                        loc_all_geoNet_stats,
                        fname_all_geoNet_stats, 
                        loc_V1A, save_stats=False,
                        fname=fname, loc=os.getcwd())

#Now create another dictionary in which the lon, lat are read from .V1A files
#which may or may not be in WGS84 coordinates. Unless something goes wrong
#this should contain all the downloaded stations
#Note fname=None crashed
event_stats2, _ = sg.statsll_from_V1A(
                     loc_V1A, save_stats=False,
                     fname=fname, loc=os.getcwd())


#Find stats that are in event_stat2 but not in event_stats1
stat_codes1 = set(event_stats1.keys())
stat_codes2 = set(event_stats2.keys())

#perform check
if  not stat_codes1.issubset(stat_codes2):
    print("Some station (lon, lat) were not read from .V1A files\n")

for stat_code in (stat_codes2 -stat_codes1):
    event_stats1[stat_code] = event_stats2[stat_code]

#write statsll for this event
with open("/".join([os.getcwd(),fname]),'w') as f:
    for stat_code in event_stats1:
        (lon, lat) = event_stats1[stat_code]
        f.write("{:<15.4f} {:^15.4f} {:^10s} \n".format(lon, lat, stat_code))

#or use the convenience function
#write_statsll(loc, fname, event_stats1)
write_stats_kml(os.getcwd(), fname_statsll.split(".")[0]+".kml", event_stats1)

final_time=time()
print("Done in {:10.1f} secs".format(final_time-init_time))

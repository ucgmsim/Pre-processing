import sys

sys.path.append("/nesi/projects/nesi00213/Pre-processing/geoNet")
import scrapeGeoNet as sg
from geoNet_file import GeoNet_File
from process import Process
import os

EVENT_SUMMARY_FILE = "20150105_175001_eventStats_2016-12-12.txt"
# EVENT_SUMMARY_FILE = "20161113_162346_V1A_2016-12-11.txt"
LOC = "/nesi/projects/nesi00213/ObservedGroundMotions/Mw5pt95_20150105_174841"
# LOC = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/event_20161113_162346"


FILE_NAMES = []
with open("/".join([LOC, EVENT_SUMMARY_FILE]), "r") as f:
    for line in f.readlines():
        line = line.strip("\n")
        FILE_NAMES.append(line[:-4])

print("\n Processing %d stations in Vol1 data ..." % len(FILE_NAMES))
# print FILE_NAMES
# sys.exit("stopping execution")

for file_name in FILE_NAMES:
    print("\n**************************")
    print("%40s" % file_name)
    print("\n**************************")
    import os

    file_loc = "/".join([LOC, "Vol1", "data"])
    station_file_name = file_name + ".V1A"
    try:
        gf = GeoNet_File(station_file_name, file_loc, vol=1)
        if gf.comp_1st.acc.size < 11.0 / gf.comp_1st.delta_t:
            print("%s has less than 11 secs of data" % file_name)
            print("skipping %s" % file_name)
            continue
        gf.comp_1st.acc -= gf.comp_1st.acc.mean()
        gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
        gf.comp_up.acc -= gf.comp_up.acc.mean()
        from scipy.signal import detrend

        gf.comp_1st.acc = detrend(gf.comp_1st.acc, type="linear")
        gf.comp_2nd.acc = detrend(gf.comp_2nd.acc, type="linear")
        gf.comp_up.acc = detrend(gf.comp_up.acc, type="linear")
        pgf = Process(gf, lowcut=0.05)
        # pgf = Process(gf, lowcut=0.05, ft=0.25)
    except Exception as e:
        print(e)
        print("%s is problematic, skipping it" % station_file_name)
        continue
        # raise

    # strip numbers infront of station file names
    # station_name = ""
    # for x in file_name:
    #    if x.isalpha(): station_name += x
    station_name = file_name.split("_")[2]

    try:
        # if station_name in stats['stat']:
        pgf.save2disk(file_loc + "/velBB/", station_name, "velBB")
        pgf.save2disk(file_loc + "/velLF/", station_name, "velLF")
        pgf.save2disk(file_loc + "/accBB/", station_name, "accBB")

        # pgf.save2disk(file_loc+"/velLF_0pt25/", station_name, 'velLF')
        # pgf.save2disk(file_loc+"/accLF_0pt25/", station_name, 'accLF')
    except Exception as e:
        print(e)
        print("Skipping this station %s\n" % file_name)
        continue
print("Done..")

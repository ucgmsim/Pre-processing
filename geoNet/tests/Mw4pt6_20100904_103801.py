import sys

sys.path.append("/nesi/projects/nesi00213/Pre-processing/geoNet")
import scrapeGeoNet as sg
from geoNet_file import GeoNet_File
from process import Process
import os


EVENT_SUMMARY_FILE = "20100904_103801.CSV"
# LOC = "/hpc/home/man56/ObservedGroundMotions/Mw4pt6_20100904_103801"
LOC = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt6_20100904_103801"
# "/".join([os.getcwd(),"tests","data"])

FILE_NAMES = sg.read_GeoNet_stat_names(LOC, EVENT_SUMMARY_FILE)

print("\n Processing Vol1 data ...")
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
    except Exception as e:
        print(e)
        raise

    # strip numbers infront of station file names
    station_name = ""
    for x in file_name:
        if x.isalpha():
            station_name += x

    try:
        import numpy as np

        #        fsouth_stats = "/hpc/home/bab70/StationInfo/southislandstations_v1.ll"
        #        fcant_stats = "/hpc/home/bab70/StationInfo/cantstations.ll"
        fcant_stats = "/nesi/projects/nesi00213/StationInfo/cantstations.ll"
        stats = np.genfromtxt(
            fcant_stats, dtype="f,f,S4", names=["lon", "lat", "stat"]
        )  # , delimiter=" ")
        if station_name in stats["stat"]:
            pgf.save2disk(file_loc + "/velBB/", station_name, "velBB")
            pgf.save2disk(file_loc + "/velLF/", station_name, "velLF")
            pgf.save2disk(file_loc + "/accBB/", station_name, "accBB")

    except Exception as e:
        print(e)
        print("Skipping this station %s\n" % file_name)
        continue

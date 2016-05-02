import scrapeGeoNet as sg
from geoNet_file import GeoNet_File
from process import Process
import os
"""
We use the valentine's day 2016 earth quake event as an example.
step one:
        Set BASE_URL
step two:
        Download and save EVENT_SUMMARY_FILE to the current directory
WARNING:
        The EVENT_SUMMARY_FILE contains some station names which should 
        not be there. Remove these stations from EVENT_SUMMARY_FILE.
        This is an error on the part of the person maintaining GeoNet
"""
BASE_URL = "ftp://ftp.geonet.org.nz/strong/processed/Proc/2016/02_Feb/2016-02-14_001343/" 
EVENT_SUMMARY_FILE = "20160214_001343.CSV"
LOC = "/".join([os.getcwd(),"tests","data"])

FILE_NAMES = sg.read_GeoNet_stat_names(LOC, EVENT_SUMMARY_FILE)
sg.make_dataPlot_dirs(num_vol=4)
print("Downloading data ...")
sg.save_GeoNet_dataPlots(FILE_NAMES, BASE_URL)
print("Finished downloading data")

print("\n Processing Vol1 data ...")
for file_name in FILE_NAMES:
    import os
    file_loc = "/".join([os.getcwd(), "Vol1", "data"])
    station_file_name = file_name + ".V1A"
    try:
        gf = GeoNet_File(station_file_name, file_loc, vol=1)
    except Exception as e:
        print(e)
        print("Skipping this station %s\n" %file_name)
        continue
    gf_processed  = Process(gf)
    
    #strip numbers infront of station file names
    station_name = ""
    for x in file_name:
        if x.isalpha(): station_name += x

    gf_processed.save2disk(file_loc, station_name)

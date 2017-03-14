import sys
sys.path.append("/nesi/projects/nesi00213/Pre-processing/geoNet")
import scrapeGeoNet as sg
from geoNet_file import GeoNet_File
from process import Process
import os

#EVENT_SUMMARY_FILE = "20161113_110256.CSV"
#EVENT_SUMMARY_FILE = "20161113_110256.txt"
#EVENT_SUMMARY_FILE = "20161113_110256_missed_stations.txt"
#EVENT_SUMMARY_FILE = "20161113_110256_19Nov2016_extra_stats.txt"
EVENT_SUMMARY_FILE= "20161113_110256_19Nov2016.txt"
#LOC = "/hpc/home/man56/ObservedGroundMotions/Mw4pt9_20110429_190804" 
LOC = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw7pt5_20161113_110256"

#"/".join([os.getcwd(),"tests","data"])

#FILE_NAMES = sg.read_GeoNet_stat_names(LOC, EVENT_SUMMARY_FILE)
FILE_NAMES = []
with open("/".join([LOC,EVENT_SUMMARY_FILE]), 'r') as f:
    for line in f.readlines():
        #print("line is %s" %line)
        #print line.strip(".V1A\n")
        #for some reason strip does not work on for example 20161113_110315_GODS_21
        #FILE_NAMES.append(line.strip(".V1A\n"))
        line = line.strip("\n")
        FILE_NAMES.append(line[:-4])

#FILE_NAMES = ["20161113_110317_KEKS_20"]
print("\n Processing %d stations in Vol1 data ..." %len(FILE_NAMES))
#print FILE_NAMES
#sys.exit("stopping execution")

for file_name in FILE_NAMES:
    print("\n**************************")
    print("%40s" %file_name)
    print("\n**************************")
    import os
    file_loc = "/".join([LOC, "Vol1", "data"])
    station_file_name = file_name + ".V1A"
    try:
        gf = GeoNet_File(station_file_name, file_loc, vol=1)
        if gf.comp_1st.acc.size < 11./gf.comp_1st.delta_t:
            print("%s has less than 11 secs of data" %file_name)
            print("skipping %s" %file_name)
            continue
        gf.comp_1st.acc -= gf.comp_1st.acc.mean()
        gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
        gf.comp_up.acc  -= gf.comp_up.acc.mean()
        from scipy.signal import detrend
        gf.comp_1st.acc = detrend(gf.comp_1st.acc, type='linear')
        gf.comp_2nd.acc = detrend(gf.comp_2nd.acc, type='linear')
        gf.comp_up.acc = detrend(gf.comp_up.acc, type='linear')       
        #pgf = Process(gf, lowcut=0.05)
        #pgf = Process(gf, lowcut=0.05, ft=0.25)
        pgf = Process(gf, lowcut=0.05, ft=0.5)
    except Exception as e:
        print(e)
        print("%s is problematic, skipping it" %station_file_name)
        #raise

    #strip numbers infront of station file names
    station_name = ""
    for x in file_name:
        if x.isalpha(): station_name += x

    try:
#        import numpy as np
#        fsouth_stats = "/hpc/home/bab70/StationInfo/southislandstations_v1.ll"
#        fcant_stats = "/hpc/home/bab70/StationInfo/cantstations.ll"
#        #fcant_stats = os.getcwd()+"/southislandstations_v1.ll"
#        fcant_stats = os.getcwd()+"/amberley.ll"
#        stats = np.genfromtxt(fcant_stats, dtype="f,f,S4",
#                                    names=['lon','lat','stat'])#, delimiter=" ")
        #if station_name in stats['stat']:
        #pgf.save2disk(file_loc+"/velBB/", station_name, 'velBB')
        #pgf.save2disk(file_loc+"/velLF/", station_name, 'velLF')
        #pgf.save2disk(file_loc+"/accBB/", station_name, 'accBB')

        #pgf.save2disk(file_loc+"/velLF_0pt25/", station_name, 'velLF')
        #pgf.save2disk(file_loc+"/accLF_0pt25/", station_name, 'accLF')

        if not os.path.exists(file_loc+"/velLF_0pt5"):
            os.mkdir(file_loc+"/velLF_0pt5")
        if not os.path.exists(file_loc+"/accLF_0pt5"):
            os.mkdir(file_loc+"/accLF_0pt5")
        pgf.save2disk(file_loc+"/velLF_0pt5/", station_name, 'velLF')
        pgf.save2disk(file_loc+"/accLF_0pt5/", station_name, 'accLF')
    except Exception as e:
        print(e)
        print("Skipping this station %s\n" %file_name)
        continue
print("Done..")

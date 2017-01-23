import sys
import os
from glob import glob
from scipy.signal import detrend
from time import time
from pyqc.geoNet import scrapeGeoNet as sg
from pyqc.geoNet.geoNet_file import GeoNet_File
from pyqc.geoNet.process import Process

init_time = time()

LOC= "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt1_20110707_063901/Vol1/data"
FILE_NAMES = []
event_stats_V1A = glob("/".join([LOC,"*.V1A"]))
FILE_NAMES = [os.path.basename(_) for _ in event_stats_V1A]

print("\n Processing %d stations in Vol1 data ..." %len(FILE_NAMES))
#print FILE_NAMES
#sys.exit("stopping execution")

for station_file_name in FILE_NAMES:
    print("\n**************************")
    print("%s" %station_file_name)
    print("\n**************************")
    try:
        gf = GeoNet_File(station_file_name, LOC, vol=1)
        if gf.comp_1st.acc.size < 11./gf.comp_1st.delta_t:
            print("%s has less than 11 secs of data" %station_file_name)
            print("skipping %s" %station_file_name)
            continue
        gf.comp_1st.acc -= gf.comp_1st.acc.mean()
        gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
        gf.comp_up.acc  -= gf.comp_up.acc.mean()

        gf.comp_1st.acc = detrend(gf.comp_1st.acc, type='linear')
        gf.comp_2nd.acc = detrend(gf.comp_2nd.acc, type='linear')
        gf.comp_up.acc = detrend(gf.comp_up.acc, type='linear')       
        pgf = Process(gf, lowcut=0.05)
        #pgf = Process(gf, lowcut=0.05, ft=0.25)
    except Exception as e:
        print(e)
        print("%s is problematic, skipping it" %station_file_name)
        continue
        #raise

    #extract stat code
    stat_code = station_file_name.split(".")[0].split("_")[2]

    try:
        pgf.save2disk(LOC+"/velBB/", stat_code, 'velBB')
        pgf.save2disk(LOC+"/velLF/", stat_code, 'velLF')
        pgf.save2disk(LOC+"/accBB/", stat_code, 'accBB')

        #pgf.save2disk(file_loc+"/velLF_0pt25/", stat_code, 'velLF')
        #pgf.save2disk(file_loc+"/accLF_0pt25/", stat_code, 'accLF')
    except Exception as e:
        print(e)
        print("Skipping this station %s\n" %station_file_name)
        continue

final_time = time()
print("Done in {:10.1f} secs".format(final_time-init_time))

import os
from glob import glob
from scipy.signal import detrend
from time import time
from geoNet.geoNet.geoNet_file import GeoNet_File
from geoNet import Process, adjust_gf_for_time_delay

init_time = time()

#enter location of ObservedGroundMotions directory
#files are placed after download
LOC=os.path.join(os.getcwd(), 'ObservedGroundMotions')

#open file with eventid
f = open('list_all')

#loop over events
for line in f:
    eventid = line.strip()
    LOC_V1A = os.path.join(LOC,eventid,'Vol1','data')
    #get list of V1A files for this event
    FILE_NAMES = []
    event_stats_V1A = glob(os.path.join(LOC_V1A,'*.V1A'))
    FILE_NAMES = [os.path.basename(_) for _ in event_stats_V1A]

    print("\n Processing %d stations in Vol1 data ..." %len(FILE_NAMES))

    #loop over each V1A file and process
    for station_file_name in FILE_NAMES:
        print("\n**************************")
        print("%s" %station_file_name)
        print("\n**************************")
        try:
            gf = GeoNet_File(station_file_name, LOC_V1A, vol=1)
            if gf.comp_1st.acc.size < 5./gf.comp_1st.delta_t:
                print("%s has less than 5 secs of data" %station_file_name)
                print("skipping %s" %station_file_name)
                continue
            
            #When appended zeroes at the beginning of the record are removed, the 
            #record might then be empty, skipp processing in such a case
            agf=adjust_gf_for_time_delay(gf)
            if agf.comp_1st.acc.size <= 10: 
                print("no elements in %s. Skipping it." %station_file_name)
                continue

            gf.comp_1st.acc -= gf.comp_1st.acc.mean()
            gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
            gf.comp_up.acc  -= gf.comp_up.acc.mean()

            gf.comp_1st.acc = detrend(gf.comp_1st.acc, type='linear')
            gf.comp_2nd.acc = detrend(gf.comp_2nd.acc, type='linear')
            gf.comp_up.acc = detrend(gf.comp_up.acc, type='linear')       
            pgf = Process(gf, lowcut=0.05, gpInt=False)
            #pgf = Process(gf, lowcut=0.05, ft=0.25)
        except Exception as e:
            print(e)
            print("%s is problematic, skipping it" %station_file_name)
            continue
            #raise

        #extract stat code
        stat_code = station_file_name.split(".")[0].split("_")[2]

        try:
            #expects that velBB etc directories already exist, created when getData.py is used
            pgf.save2disk(LOC_V1A+"/velBB/", stat_code, 'velBB')
            pgf.save2disk(LOC_V1A+"/velLF/", stat_code, 'velLF')
            pgf.save2disk(LOC_V1A+"/accBB/", stat_code, 'accBB')

        except Exception as e:
            print(e)
            print("Skipping this station %s\n" %station_file_name)
            continue

final_time = time()
print("Done in {:10.1f} secs".format(final_time-init_time))

import os
from time import time
from geoNet import scrapeGeoNet as sg
from geoNet.gen_stats_kml import write_stats_kml

init_time=time()

#only the BASE_RUL needs to be changed in this example
BASE_URL="ftp://ftp.geonet.org.nz/strong/processed/Proc/2013/05_May/2013-05-04_145921/Vol1/data/"

sg.make_dataPlot_dirs_v2(loc=".",num_vol=1)
print("Downloading data ...")
LOC = "/".join([os.getcwd(),"Vol1"])

std_out, std_err = sg.get_geoNet_data(loc=LOC,
                   geoNet_dir="data",
                   geoNet_url=BASE_URL,
                   wget_options="-N"
                   )

print("Finished downloading data")
with open("std_out.txt", 'w') as f:
    f.writelines(std_out)

with open("std_err.txt", 'w') as f:
    f.writelines(std_err)



event_stats_dict, fname_statsll=sg.get_events_stations(
            fname_all_geoNet_stats="all_geoNet_stats+2016-12-20.ll",
            loc_all_geoNet_stats=os.getcwd(),
            loc_Vol1="/".join([os.getcwd(), "Vol1"]),
            save_stats=True,
            fname=None,
            loc=os.getcwd())


#remove stations for which don't know the lon, lat
for stat_code in event_stats_dict.keys():
    if event_stats_dict[stat_code][0] is None:
        event_stats_dict.pop(stat_code)

write_stats_kml(os.getcwd(), fname_statsll.split(".")[0]+".kml", event_stats_dict)

final_time=time()
print("Done in {:10.1f} secs".format(final_time-init_time))


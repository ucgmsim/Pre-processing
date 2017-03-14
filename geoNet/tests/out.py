import os
#import scrapeGeoNet as sg
from pyqc.geoNet import scrapeGeoNet as sg

BASE_URL = "ftp://ftp.geonet.org.nz/strong/processed/Proc/2015/01_Jan/2015-01-05_174841/Vol1/data"

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



event_stats=sg.get_events_stations(
            fname_all_geoNet_stats="all_geoNet_stats.ll",
            loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo",
            loc_Vol1="/".join([os.getcwd(), "Vol1"]),
            save_stats=True,
            fname=None,
            loc=os.getcwd())

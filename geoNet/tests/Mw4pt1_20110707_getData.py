import os
from time import time
from pyqc.geoNet import scrapeGeoNet as sg

# import scrapeGeoNet as sg

init_time = time()
# BASE_URL = "ftp://ftp.geonet.org.nz/strong/processed/Proc/2015/01_Jan/2015-01-05_174841/Vol1/data"
BASE_URL = "ftp://ftp.geonet.org.nz/strong/processed/Proc/2011/07_Jul/2011-07-07_063901/Vol1/data"

# sg.make_dataPlot_dirs_v2(loc=".",num_vol=1)
# print("Downloading data ...")
# LOC = "/".join([os.getcwd(),"Vol1"])
#
# std_out, std_err = sg.get_geoNet_data(loc=LOC,
#                   geoNet_dir="data",
#                   geoNet_url=BASE_URL,
#                   wget_options="-N"
#                   )
#
# print("Finished downloading data")
# with open("std_out.txt", 'w') as f:
#    f.writelines(std_out)
#
# with open("std_err.txt", 'w') as f:
#    f.writelines(std_err)
#


event_stats = sg.get_events_stations(
    fname_all_geoNet_stats="all_geoNet_stats+2016-12-20.ll",
    loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo",
    loc_Vol1="/".join([os.getcwd(), "Vol1"]),
    save_stats=True,
    fname=None,
    loc=os.getcwd(),
)

final_time = time()
print("Done in {:10.1f} secs".format(final_time - init_time))

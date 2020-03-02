import os
from time import time
from geoNet.geoNet import scrapeGeoNet as sg

init_time=time()

years = ['2019']

for i,x in enumerate(years):
    BASE_URL='ftp://ftp.geonet.org.nz/strong/processed/' + x
    
    std_out, std_err = sg.get_geoNet_data(loc=os.getcwd(),
                       geoNet_dir=os.getcwd(),
                       geoNet_url=BASE_URL,
                       wget_options="-N -t inf -r -l10"
                       )
    
    print("Finished downloading data")
    with open("std_out.txt", 'w') as f:
        f.writelines(std_out)
    
    with open("std_err.txt", 'w') as f:
        f.writelines(std_err)

final_time=time()
print("Done in {:10.1f} secs".format(final_time-init_time))

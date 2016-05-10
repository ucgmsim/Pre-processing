import numpy as np
import os
from os import system as sys
import csv
import urllib2


def read_GeoNet_stat_names(loc, id_summary_file):
    f = open("/".join([loc,id_summary_file]),'r')
    fcsv = csv.DictReader(f,delimiter=",")#fcsv is an iterator
    FILE_NAMES = [ _['Accelerogram ID'] for _ in fcsv ] 
    f.close()

    return FILE_NAMES

#create directory structure for earthquake event
def make_dataPlot_dirs(num_vol=4):
    for i in range(num_vol):
        dirName = "Vol"+str(i+1)
        os.mkdir(dirName)
        os.mkdir("/".join([dirName, "data"]))
        os.mkdir("/".join([dirName, "plots"]))
    
    return
#save file names to be downloaded in each directory and subdirectory
cwd = os.getcwd()
def save_GeoNet_dataPlots(FILE_NAMES, BASE_URL, num_vol=4):
    for i in range(1,num_vol+1):
        if i in [3,4]: continue
        vol = "Vol"+str(i)
        data_ext = "V%dA" % i

        for station in FILE_NAMES:
            url = "/".join([BASE_URL,vol,"data",station + "." + data_ext])
            stat_data = None
            try:
                stat_data = urllib2.urlopen(url)
            except Exception as  e:
                print("Remove station field below from .CSV file")
                print(e)
                continue

            with open(".".join([station,data_ext]),"w") as f:
                f.writelines(stat_data.readlines())
        
            sys("mv " + ".".join([station,data_ext]) + " "
                      +  "/".join([vol,"data"]))
    
    return


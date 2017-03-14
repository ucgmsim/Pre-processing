import numpy as np
import os
from os import system as sys
import csv
import urllib2
from glob import glob
import datetime
import subprocess as sp

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

def mkdir_p(dirName):
    if not os.path.exists(dirName):
        os.mkdir(dirName)

    return

def make_dataPlot_dirs_v2(loc=os.getcwd(), num_vol=4):
    """
    os.makedirs is also an option, but if the user makes a typo then all the 
    directories will be created. Better to use mkdir
    loc:
        parent directory where Vol1 will be created
    """
    
    for i in range(num_vol):
        dirName = "Vol"+str(i+1)
        dirName = "/".join([loc,dirName])
        mkdir_p(dirName)
        mkdir_p("/".join([dirName, "data"]))
        mkdir_p("/".join([dirName, "plots"]))
        mkdir_p("/".join([dirName, "data", "accBB"]))
        mkdir_p("/".join([dirName, "data", "velBB"]))
        mkdir_p("/".join([dirName, "data", "dispBB"]))
        mkdir_p("/".join([dirName, "data", "accLF"]))
        mkdir_p("/".join([dirName, "data", "velLF"]))
        mkdir_p("/".join([dirName, "data", "dispLF"]))
    
    return




def get_geoNet_data(loc=None, geoNet_dir='data', geoNet_url=None, wget_options='-N'):
    """
    Wrapper function that calls wget.

    loc:
        directory where the geoNet data will be saved
    outDir: 
        Directory where acceleration seismograms are placed
    geoNet_dir:
        name of geoNet directory tobe downloaded
    """
    assert loc is not None, "Specify location where geoNet data will be placed"
    assert geoNet_url is not None, "Specify geoNet url to download"

    cwd = os.getcwd()
    DirData=loc
    os.chdir(DirData)

    Prog = "wget"
    wget_default_options = " -r -np -nd -nH "
    wget_options += wget_default_options
    wget_input = "%s %s %s -P %s" %(Prog, wget_options, geoNet_url, geoNet_dir)
    wget = sp.Popen(wget_input, 
                  stdin  = sp.PIPE,
                  stdout = sp.PIPE,
                  stderr = sp.PIPE,
                  shell=True) #stderr = sp.STDOUT)

    stdout, stderr = wget.communicate()

    os.chdir(cwd)

    return stdout, stderr


def read_statsll(loc, name):
    """
    fname:
        .ll format station file that contains lon, lat, stat_code
    returns a dictionary
    """
    #lon, lat, stat_code = np.genfromtxt(fname, dtype="%10.4f, %10.4f, %5s", unpack=True)
    #stats_dict = {"lon":lon, "lat":lat, "stat":stat_code}

    stats_dict = {}
    with open("/".join([loc, name]), 'r') as f:
        for line in f:
            lon, lat, stat_code = line.split()
            lon = float(lon)
            lat = float(lat)
            stat_code = stat_code.strip() #remove space and newline character
            stats_dict[stat_code] = (lon, lat)

    return stats_dict

def statcode_lonlat_from_V1A(loc,fname):
    """
    loc: location of where geoNet data file (.V1A) is located
    fname: name of .V1A file
    Note: See http://info.geonet.org.nz/display/appdata/Accelerogram+Data+Filenames+and+Formats
        for geoNet file formats
    """
    with open ("/".join([loc,fname]), 'r') as f:
        for line_num, line in enumerate(f):
            if line_num == 1: line2 = line
            if line_num == 18:
                line19 = line
                break

    #assume line 2 is of the form:
    #Site CTZ       43 44 14S  176 37 03W
    stat_code = line2.split()[1]
    line19 = line19.split()
    lat = -(float(line19[0]) + float(line19[1])/60 + float(line19[2])/3600)
    lon = float(line19[3]) + float(line19[4])/60 + float(line19[5])/3600

    return {stat_code: (lon, lat)} #{stat_code: (round(lon,4), round(lat,4))}

def parse_V1A_fname(fname):
    """
    Note: See http://info.geonet.org.nz/display/appdata/Accelerogram+Data+Filenames+and+Formats
    The file name can take three forms:
    XNNSITEJJ, where X is the instrument code1, NN is the two-digit year 
    of acquisition, SITE is the 4 character site code, and JJ is the site's
    unique accelerogram number for that year.
    YYYYMMDD_HHMM_SITE, where YYYY/MM/DD HH:MM is the earthquake origin time 
    (UTC) and SITE is the 3 or 4 character site code.
    YYYYMMDD_HHMM_SITE_XX, where YYYY/MM/DD HH:MM is the instrument trigger time
    (UTC), SITE is the 3 or 4 character site code, and XX is the sensor location
    code. The latter is most appropriate for arrays of borehole sensors.

    fname: name of geoNet data file, must not include path
    """
    form1 = "XNNSITEJJ"
    form2 = "YYYYMMDD_HHMM_SITE"
    form3 = "YYYYMMDD_HHMM_SITE_XX"

    #Remove .V1A, .V2A extension
    fname = fname.split(".")[0]
    YYYYMMDD = ""
    HHMM = ""
    if fname.count("_") == form1.count("_"):
        stat_code = fname[3:7]
    elif fname.count("_") == form2.count("_"):
        YYYYMMDD, HHMM, stat_code = fname.split("_")
    elif fname.count("_") == form3.count("_"):
        YYYYMMDD, HHMM, stat_code, XX = fname.split("_")
    else:
        raise Exception("{:s} is unknow file name format for .V1A files\n".format(fname))

    return YYYYMMDD, HHMM, stat_code

def statsll_from_V1A(loc_V1A, save_stats=False,
                     fname=None, loc=os.getcwd()):
    """
    loc_V1A: location where .V1A geoNet data fils are located
    save_stats: If save_stats==True, save  fname.ll at loc. Note that 
    fname and loc are not used if save_stats==False
    Note:
        Here (lon, lat) are extracted from .V1A file
    return:
        (1) A dictionary containing station code, lon and lat in the format of 
        read_statsll. 
        (2) station file name

    """
    event_stats_V1A = glob("/".join([loc_V1A,"*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    event_stats = {}
    for V1A_file in event_stats_V1A:
        try:
            event_stats.update(
            statcode_lonlat_from_V1A(loc_V1A, V1A_file)
            )
        except Exception as e:
            print(e)
            print("Stations list will not include {:s}\n".format(V1A_file))
            continue

    #year = ""
    #event_id = ""
    for V1A_file in event_stats_V1A:
        try:
            #year, event_id, stat_code, V1A = V1A_file.split(".")[0].split("_")
            #split_file_name = V1A_file.split(".")[0].split("_")
            #year, event_id, stat_code = split_file_name[0:3]
            YYYYMMDD, HHMM, stat_code = parse_V1A_fname(V1A_file)
        except Exception as e:
            print(e)
            continue
        else:
            if YYYYMMDD is "":
                continue
            else:
                break

    if fname is None:
        fname = "_".join([YYYYMMDD, "eventStats", str(datetime.date.today())])
        fname+=".ll"

    if save_stats is True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        with open("/".join([loc,fname]),'w') as f:
            for stat_code, (lon, lat) in event_stats.items():
                f.write("{:^10.4f} {:^10.4f} {:^10s} \n".format(lon, lat, stat_code))

    return event_stats, fname


def statsll_from_V1A_v2(loc_all_geoNet_stats, fname_all_geoNet_stats, 
                        loc_V1A, save_stats=False,
                        fname=None, loc=os.getcwd()):
    """
    Note: This functions assigns (lon, lat) values after reading them from 
    loc_all_geoNet_stats. If save_stats is True then fname should be suppllied.
    loc_all_geoNet_stats must be in WGS84 coordinates as in 
    https://magma.geonet.org.nz/delta/app
    return:
        (1) A dictionary containing station code, lon and lat in the format of 
        read_statsll. 
    """

    all_geoNet_stats = read_statsll(loc_all_geoNet_stats, fname_all_geoNet_stats)

    event_stats_V1A = glob("/".join([loc_V1A,"*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    event_stats = {}
    stats_not_in_all_geoNet_stats=[]
    for V1A_file in event_stats_V1A:
        YYYYMMDD, HHMM, stat_code = parse_V1A_fname(V1A_file)
        if all_geoNet_stats.has_key(stat_code):
            event_stats[stat_code] = all_geoNet_stats[stat_code]
        else:
            stats_not_in_all_geoNet_stats.append(stat_code)

    if len(stats_not_in_all_geoNet_stats) > 0:
        print("(lon, lat) for the following stations were not in {:s}:\n".format(fname_all_geoNet_stats))
        print(stats_not_in_all_geoNet_stats)

    if fname is None:
        fname = "_".join([YYYYMMDD, "eventStats", str(datetime.date.today())])
        fname+=".ll"

    if save_stats is True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        with open("/".join([loc,fname]),'w') as f:
            for stat_code, (lon, lat) in event_stats.items():
                f.write("{:^10.4f} {:^10.4f} {:^10s} \n".format(lon, lat, stat_code))

    return event_stats, fname


def get_events_stations(fname_all_geoNet_stats=None,
                        loc_all_geoNet_stats=None,
                        loc_Vol1="/".join([os.getcwd(), "Vol1"]),
                        save_stats=False, fname=None, loc=os.getcwd()):
    """
    This function should no longer be used
    """

    all_geoNet_stats = read_statsll(loc_all_geoNet_stats, fname_all_geoNet_stats)

    event_stats_V1A = glob("/".join([loc_Vol1,"data","*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    
    event_stats= {}
    for V1A_file in event_stats_V1A:
        #year, event_id, stat_code, V1A = V1A_file.split(".")[0].split("_")
        split_file_name = V1A_file.split(".")[0].split("_")
        year, event_id, stat_code = split_file_name[0:3]
        event_stats[stat_code] = (None, None)
        if all_geoNet_stats.has_key(stat_code):
            event_stats[stat_code] = all_geoNet_stats[stat_code]

    if save_stats == True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        if fname is None:
            fname = "_".join([year, event_id, "eventStats", str(datetime.date.today())])
            fname+=".ll"
            with open("/".join([loc,fname]),'w') as f:
                for key, value in event_stats.items():
                    if value[0] is None:
                        print("{:10s} not found in all_geoNet_stats, add this manually to event_stats.ll".format(key))
                    else:
                        line="{:10.4f} {:10.4f} {:10s}".format(value[0], value[1], key)
                    f.write(line+"\n")

    return event_stats, fname







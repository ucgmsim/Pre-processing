"""
Reads event_info.txt and performs task for real-time simulation drills.
The functions here are named after and perform the same tasks as those 
separate python scripts. This python script is specifically written for use 
in real-time simulations drills. Each function performs a specific and indepedent 
task. For individual needs and customized interaction with geoNet package use 
the separate python scripts instead.
"""

import numpy as np
import os
from time import time
from glob import glob
import scipy
import sys
import subprocess as sp
#geoNet imports
from geoNet import scrapeGeoNet as sg
from geoNet import utils, putils
from geoNet.gen_stats_kml import write_stats_kml
from geoNet.geoNet_file import GeoNet_File
from geoNet.process import Process, adjust_gf_for_time_delay
import create_1d_profile


def getData(loc, BASE_URL):
    """
    loc:
        .V1A data will be downloaded to loc/Vol1
    """
    init_time=time()
    print("Downloading data ...")

    sg.make_dataPlot_dirs_v2(loc=loc,num_vol=1)

    std_out, std_err = sg.get_geoNet_data(
                       loc="/".join([loc, "Vol1"]),
                       geoNet_dir="data",
                       geoNet_url=BASE_URL,
                       wget_options="-N"
                       )

    with open("getData_stdout.txt", 'w') as f:
        f.writelines(std_out)

    with open("getData_stderr.txt", 'w') as f:
        f.writelines(std_err)


    final_time=time()
    print("Finished downloading data in {:.1f} secs\n".format(final_time - init_time))

    return


def event_statsll(fname, loc,
                  loc_all_geoNet_stats, fname_all_geoNet_stats,
                  loc_V1A):
    """
    fname:
        stations.ll file name
    loc:
        stations.ll saved at loc
    loc_all_geoNet_stats, fname_all_geoNet_stats:
        location and name of master file that contains (lon, lat, stat_code)
        for all SMSs.
    loc_V1A:
        where .V1A data files are saved in getData.
    """
    init_time=time()
    print("Creating stations list ...")
    
    std_out = open("event_statsll_stdout.txt", 'w')
    #Firs get statsll dictionary by assigning lon, lat from the saved list known 
    #to be in WGS84 coordinates
    (
    event_stats1,
    fname_statsll
    ) = sg.statsll_from_V1A_v2(
                            loc_all_geoNet_stats,
                            fname_all_geoNet_stats, 
                            loc_V1A, save_stats=False,
                            fname=fname, loc=loc)

    #Now create another dictionary in which the lon, lat are read from .V1A files
    #which may or may not be in WGS84 coordinates. Unless something goes wrong
    #this should contain all the downloaded stations
    #Note fname=None crashed
    event_stats2, _ = sg.statsll_from_V1A(
                         loc_V1A, save_stats=False,
                         fname=None, loc=loc)


    #Find stats that are in event_stat2 but not in event_stats1
    stat_codes1 = set(event_stats1.keys())
    stat_codes2 = set(event_stats2.keys())

    #perform check
    if  not stat_codes1.issubset(stat_codes2):
        std_out.write("Some station (lon, lat) were not read from .V1A files\n")

    for stat_code in (stat_codes2 -stat_codes1):
        event_stats1[stat_code] = event_stats2[stat_code]
    #write statsll for this event
    fpath=os.path.join(loc,fname)
    if not os.path.isdir(loc):
        os.makedirs(loc)
        
    with open(fpath,'w') as f:
        for stat_code in event_stats1:
            (lon, lat) = event_stats1[stat_code]
            f.write("{:<15.4f} {:^15.4f} {:^10s} \n".format(lon, lat, stat_code))

    #write kml file to view with google earth/map
    write_stats_kml(loc, fname_statsll.split(".")[0]+".kml", event_stats1)

    std_out.close()
    final_time=time()
    print("Done in {:.1f} secs\n".format(final_time - init_time))

    return


def event_statsVs30(fname_statsll, loc=os.getcwd(),
                  dirVs30_prog="/home/ahsan/Documents/Vs30-mapping",
                  Vs30_prog="extractVs30.R"):
    """
    fname_statsll:
        stations.ll file name, pass the file name and full path combined
    loc:
        where .Vs30 files are saved 
    """
    init_time=time()
    print("Creating Vs30 files  ...")

    cwd = os.getcwd()
    dirProg = dirVs30_prog
    os.chdir(dirProg)


    fname_Vs30_out = "/".join([loc,"Vs30.out"])

    fname_Vs30_in="/".join([loc, "statsll.temp"])
    with open(fname_statsll, 'r') as f, open(fname_Vs30_in, 'w') as g:
        g.write("Site.Longitude  Site.Latitude  Site.Code\n")
        for line in f:
            g.write(line)

    prog = sp.Popen(["Rscript", Vs30_prog, "KRIGE_NZGD00_allNZ.Rdata", fname_Vs30_in,
                    fname_Vs30_out], stdout=sp.PIPE,
                    stderr=sp.PIPE, shell=False)

    std_out, std_err = prog.communicate()
    with open("/".join([loc, "Vs30_stdout.txt"]), 'w') as f:
        f.writelines(std_out)

    with open("/".join([loc, "Vs30_stderr.txt"]), 'w') as f:
        f.writelines(std_err)

    os.remove(fname_Vs30_in)

    fname_Vs30 = fname_statsll.replace(".ll", ".vs30")
    fname_Vs30Ref = fname_statsll.replace(".ll", ".vs30ref")
    with open(fname_Vs30_out, 'r') as fin, open(fname_Vs30, 'w') as fVs30, open(fname_Vs30Ref, 'w') as fVs30_ref:
        header = fin.readline()
        fVs30.write("%used {:s}\n".format("/".join([dirVs30_prog, Vs30_prog])))
        fVs30_ref.write("%used {:s}\n".format("/".join([dirVs30_prog, Vs30_prog])))
        for line in fin:
            stat_code, lon, lat, Vs30 = line.split()
            #remove \n
            Vs30 = Vs30.strip()
            stat_code = stat_code.strip("\"")
            if Vs30 == "NA":
                Vs30 = 500

            Vs30 = float(Vs30)

            fVs30.write("{:<10s} {:^10f}\n".format(stat_code, Vs30))
            fVs30_ref.write("{:<10s} {:^10d}\n".format(stat_code, 500))
            #lon=float(lon)
            #lat=float(lat)
            #fout.write("{:<15.4f} {:^15.4f} {:^10d}\n".format(lon, lat, Vs30))       
            #fout.write("{:<15.4f} {:^15.4f} {:^10s}\n".format(lon, lat, stat_code))
            #fout.write("{:<15.4f} {:^15.4f} {:^10d} {:^10s}\n".format(lon, lat, Vs30, stat_code))

    os.chdir(cwd)
    final_time=time()
    print("Done Vs30 calcs in {:.1f} secs\n".format(final_time - init_time))

    return


def processData(LOC):
    """
    LOC:
        where .V1A data files are saved in getData.
    """
    init_time=time()
    print("Processing .V1A files ...")

    std_err = open("processData_stderr.txt", 'w')
    std_out = open("processData_stdout.txt", 'w')
    sys.stdout = std_out


    FILE_NAMES = []
    event_stats_V1A = glob("/".join([LOC,"*.V1A"]))
    FILE_NAMES = [os.path.basename(_) for _ in event_stats_V1A]

    std_out.write("\n Processing %d stations in Vol1 data ...\n" %len(FILE_NAMES))

    for station_file_name in FILE_NAMES:
        std_out.write("\n**************************\n")
        std_out.write("%s" %station_file_name)
        std_out.write("\n**************************\n")
        try:
            gf = GeoNet_File(station_file_name, LOC, vol=1)
            if gf.comp_1st.acc.size < 5./gf.comp_1st.delta_t:
                std_out.write("%s has less than 5 secs of data\n" %station_file_name)
                std_out.write("skipping %s\n" %station_file_name)
                continue
            
            #When appended zeroes at the beginning of the record are removed, the 
            #record might then be empty, skipp processing in such a case. 
            #The number 10 is rather arbirary, 1 should suffice
            agf=adjust_gf_for_time_delay(gf)
            if agf.comp_1st.acc.size <= 10: 
                std_out.write("no elements in %s. Skipping it.\n" %station_file_name)
                continue

            gf.comp_1st.acc = scipy.signal.detrend(gf.comp_1st.acc, type='linear')
            gf.comp_2nd.acc = scipy.signal.detrend(gf.comp_2nd.acc, type='linear')
            gf.comp_up.acc = scipy.signal.detrend(gf.comp_up.acc, type='linear')       

            gf.comp_1st.acc -= gf.comp_1st.acc.mean()
            gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
            gf.comp_up.acc  -= gf.comp_up.acc.mean()


            pgf = Process(gf, lowcut=0.05, gpInt=False)
            #pgf = Process(gf, lowcut=0.05, ft=0.25)
        except Exception as e:
            std_err.write(str(e))
            std_err.write("%s is problematic, skipping it.\n" %station_file_name)
            continue
            #raise

        #extract stat code
        #stat_code = station_file_name.split(".")[0].split("_")[2]
        #This could be placed in a try, except block
        YYYYMMDD, HHMM, stat_code = utils.parse_V1A_fname(station_file_name)

        try:
            #expects that velBB etc directories already exist, created when getData.py is used
            pgf.save2disk(LOC+"/velBB/", stat_code, 'velBB',comment="observed")
            pgf.save2disk(LOC+"/velLF/", stat_code, 'velLF',comment="observed")
            pgf.save2disk(LOC+"/accBB/", stat_code, 'accBB',comment="observed")

        except Exception as e:
            std_err.write(str(e))
            std_err.write("Skipping this station %s\n" %station_file_name)
            continue

    std_out.close()
    std_err.close()

    sys.stdout=sys.__stdout__
    final_time=time()
    print("Done processing in {:.1f} secs\n".format(final_time - init_time))

    return

def plot_accvel(parent_dir_loc, plot_dir_accBB, plot_dir_velBB,
                loc_statsll, fname_statsll):
    """
    Plots all the SMS data in a given directory by reading all *.000 files and 
    sorting the SMSs in order of increasing PGV
    """
    init_time=time()
    print("Creating acc and vel plots for all Observed SMSs ...")

    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
    from geoNet.utils import get_sorted_stats_code,  read_statsll, get_processed_stats_list

    init_time = time()
    loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB
    loc_vel = loc_velBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)

    #stations are sorted according to PGV
    sorted_stats_code = get_sorted_stats_code(loc_velBB,stats_dict,comp='geom')
    pdf_acc = PdfPages('plots_acc.pdf')
    pdf_vel = PdfPages('plots_vel.pdf')
    for stat in sorted_stats_code:
        stat_acc, fig_acc, ax_acc = get_stat_acc_plot(loc_acc, stat["name"])
        pdf_acc.savefig()

        stat_acc, fig_vel, ax_vel = get_stat_vel_plot(loc_vel, stat["name"])
        pdf_vel.savefig()
        plt.close('all')

    pdf_acc.close()
    pdf_vel.close()

    final_time = time()
    print("Done plotting acc, vel in {:.1f} secs.\n".format(final_time - init_time))
    return


def plot_psa(parent_dir_loc, plot_dir_accBB, plot_dir_velBB,
             loc_statsll, fname_statsll):
    """
    """
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from geoNet.putils import get_stat_PSA_plot
    from geoNet.utils import get_sorted_stats_code,  read_statsll, get_processed_stats_list

    init_time=time()
    print("Creating pSA plots for all Observed SMSs ...")

    loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)
    #stations are sorted according to PGV
    sorted_stats_code = get_sorted_stats_code(loc_velBB,stats_dict)
    pdf_psa = PdfPages('plots_psa.pdf')
    for stat in sorted_stats_code:
        stat, PSA, fig_psa, ax = get_stat_PSA_plot(loc_acc, stat["name"])
        pdf_psa.savefig()

        plt.close('all')

    pdf_psa.close()

    final_time = time()
    print("Done pSA plots in {:.1f} secs.\n".format(final_time - init_time))

    return


def IMsOnMap(parent_dir_loc, plot_dir_accBB, plot_dir_velBB,
             loc_statsll, fname_statsll):
    """
    Writes PGV, PGA, pSA intensity measures to files used in GMT plotting scripts
    """

    from geoNet.utils import (get_sorted_stats_code, read_statsll, 
    get_event_data, get_event_PSA, get_extremum, get_processed_stats_list)

    init_time=time()
    print("Writting IMs for GMT plotting ...")

    loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB
    loc_vel = loc_velBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)

    #stations are sorted according to PGV
    #Altough not required helps identify maximum and minimum IM interval
    sorted_stats_code = get_sorted_stats_code(loc_velBB, stats_dict, comp='geom')
    sorted_stats_code = [_["name"] for _ in sorted_stats_code]

    #write PGVs
    event_velBB = get_event_data(loc_velBB, sorted_stats_code)
    with open("event_obs_PGVs_lin.txt", 'w') as f, open("event_obs_PGVs_log.txt", 'w') as g:
        f.write("Observed PGVs\n")
        f.write("units of cm/s\n")
        f.write("hot:invert 0.2\n")
        f.write("0.008471  20.5 0.005 0.25\n")
        f.write("1\n")
        f.write("{:s}\n".format(fname_statsll[0:17]))
        g.write("Observed PGVs\n")
        g.write("units of cm/s\n")
        g.write("hot:invert 0.2\n")
        g.write("0.008471  20.5 0.005 0.25\n")
        g.write("1\n")
        g.write("{:s}\n".format(fname_statsll[0:17]))
        for stat_data in event_velBB:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            ext_vel_000, _ = get_extremum(stat_data['000'])
            ext_vel_090, _ = get_extremum(stat_data['090'])
            max_vel_geom = np.sqrt(np.abs(ext_vel_000*ext_vel_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line+= " {:^15.6f}".format(max_vel_geom)
            f.write(line+"\n")
            g.write("{:10.4f} {:10.4f} {:^15.6f}\n".format(lon, lat, np.log(max_vel_geom)))

    #write PGAs
    event_accBB = get_event_data(loc_accBB, sorted_stats_code)
    #with open("event_obs_PGAs.txt", 'w') as f:
    with open("event_obs_PGAs_lin.txt", 'w') as f, open("event_obs_PGAs_log.txt", 'w') as g:
        f.write("Observed PGAs\n")
        f.write("units of g\n")
        f.write("hot:invert 0.2\n")
        f.write("0.008471  20.5 0.005 0.25\n")
        f.write("1\n")
        f.write("{:s}\n".format(fname_statsll[0:17]))
        g.write("Observed PGAs\n")
        g.write("units of g\n")
        g.write("hot:invert 0.2\n")
        g.write("0.008471  20.5 0.005 0.25\n")
        g.write("1\n")
        g.write("{:s}\n".format(fname_statsll[0:17]))
        for stat_data in event_accBB:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            ext_acc_000, _ = get_extremum(stat_data['000'])
            ext_acc_090, _ = get_extremum(stat_data['090'])
            max_acc_geom = np.sqrt(np.abs(ext_acc_000*ext_acc_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line+= " {:^15.6f}".format(max_acc_geom)
            f.write(line+"\n")
            g.write("{:10.4f} {:10.4f} {:^15.6f}\n".format(lon, lat, np.log(max_acc_geom)))
    #sys.exit("getting out")


    #need new iterator object
    #write PSAs
    period=np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])
    event_PSA   = get_event_PSA(get_event_data(loc_accBB, sorted_stats_code),
                                period, xi=0.05, m=1., gamma=0.5, beta=0.25)

    #with open("event_obs_PSAs.txt", 'w') as f:
    with open("event_obs_PSAs_lin.txt", 'w') as f, open("event_obs_PSAs_log.txt", 'w') as g:
        f.write("Observed PSAs\n")
        f.write("units of g\n")
        f.write("hot:invert 0.2\n")
        f.write("0.008471  20.5 0.005 0.25\n")
        f.write("{:d}\n".format(len(period)))
        f.write((len(period)*" {:5.1f} ").format(*period))
        f.write("\n")
        g.write("Observed PSAs\n")
        g.write("units of g\n")
        g.write("hot:invert 0.2\n")
        g.write("0.008471  20.5 0.005 0.25\n")
        g.write("{:d}\n".format(len(period)))
        g.write((len(period)*" {:5.1f} ").format(*period))
        g.write("\n")

        for stat_data in event_PSA:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            pSA_geom = stat_data['geom']
            line = "{:10.4f} {:10.4f}".format(lon, lat)
            logline=line+ (len(period)*" {:^15.6f} ").format(*np.log(pSA_geom))
            line+= (len(period)*" {:^15.6f} ").format(*pSA_geom)
            f.write(line+"\n")
            g.write(logline+"\n")


    final_time = time()
    print("Done in {:.1f} secs.\n".format(final_time - init_time))
    return

#def keyValueFromTxt(fname):
#    """
#    Parses file that has the form key=value and returns a dictionary
#    """
#    keyValue = dict()
#    fname = os.path.abspath(fname)
#    print("Reading input from {:s}\n".format(fname))
#    with open(fname, 'r') as f:
#        for line in f:
#            #remove white spaces
#            if line.startswith("#") or line.startswith("%"):
#                continue
#            if line in ["\n", "\r", "\rn"]:
#                continue
#            line = line.strip()
#            line = line.replace(" ", "")
#            line = line.replace("\"", "")
#            key, value = line.split("=")
#            keyValue[key] = value
#
#    return keyValue

def run(arg):
    """
    main function that runs realtime.py
    """
    keyValue = utils.keyValueFromTxt(arg) 

    loc = keyValue['loc'] 
    BASE_URL = keyValue['BASE_URL']
    loc_statsll = keyValue['loc_statsll']
    fname_statsll = keyValue['fname_statsll']
    loc_all_geoNet_stats = keyValue['loc_all_geoNet_stats']
    fname_all_geoNet_stats = keyValue['fname_all_geoNet_stats']
    loc_V1A = keyValue['loc_V1A']
    obs_velDir = keyValue['obs_velDir']
    obs_accDir = keyValue['obs_accDir']

    site_specific= True


    if keyValue.has_key('dirVs30_prog'):
        dirVs30_prog=keyValue['dirVs30_prog']
    else:
        dirVs30_prog="/home/ahsan/Documents/Vs30-mapping"

    if keyValue.has_key('Vs30_prog'):
        Vs30_prog=keyValue['Vs30_prog']
    else:
        Vs30_prog="extractVs30.R"

    print dirVs30_prog


    getData(loc, BASE_URL)

    event_statsll(fname_statsll, loc_statsll,
                  loc_all_geoNet_stats, fname_all_geoNet_stats,
                  loc_V1A)


    event_statsVs30("/".join([loc_statsll, fname_statsll]),
                    loc=loc_statsll,
                    dirVs30_prog=dirVs30_prog,
                    Vs30_prog=Vs30_prog)

    if site_specific:
        create_1d_profile.event_generate_multiple_profiles(os.path.join(loc_statsll,fname_statsll))

    processData(loc_V1A)

    plot_accvel(loc_V1A, obs_accDir, obs_velDir,
                loc_statsll, fname_statsll)

    plot_psa(loc_V1A, obs_accDir, obs_velDir,
             loc_statsll, fname_statsll)



    IMsOnMap(loc_V1A, obs_accDir, obs_velDir,
             loc_statsll, fname_statsll)

    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Real-time simulation event_info')
    parser.add_argument('-f','--fname', help='name of input file is required', required=True)
    args = vars(parser.parse_args())
    run(args['fname'])


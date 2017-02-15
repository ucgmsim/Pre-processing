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
#geoNet imports
from geoNet import scrapeGeoNet as sg
from geoNet import utils, putils
from geoNet.gen_stats_kml import write_stats_kml
from geoNet.geoNet_file import GeoNet_File
from geoNet.process import Process, adjust_gf_for_time_delay

def getData(loc, BASE_URL):
    """
    loc:
        .V1A data will be downloaded to loc/Vol1
    """
    init_time=time()
    print("Downloading data ...")
    #only the BASE_RUL needs to be changed in this example
    #BASE_URL="ftp://ftp.geonet.org.nz/strong/processed/Proc/2017/02_Feb/2017-02-01_102129/Vol1/data/"

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
    #fname="event_stats.ll"
    ##loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo"
    #loc_all_geoNet_stats="."
    #fname_all_geoNet_stats="all_geoNet_stats+2016-12-20.ll"
    #loc_V1A="/".join([os.getcwd(), "Vol1", "data"])
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
                         fname=fname, loc=loc)


    #Find stats that are in event_stat2 but not in event_stats1
    stat_codes1 = set(event_stats1.keys())
    stat_codes2 = set(event_stats2.keys())

    #perform check
    if  not stat_codes1.issubset(stat_codes2):
        std_out.write("Some station (lon, lat) were not read from .V1A files\n")

    for stat_code in (stat_codes2 -stat_codes1):
        event_stats1[stat_code] = event_stats2[stat_code]
    #write statsll for this event
    with open("/".join([loc, fname]),'w') as f:
        for stat_code in event_stats1:
            (lon, lat) = event_stats1[stat_code]
            f.write("{:<15.4f} {:^15.4f} {:^10s} \n".format(lon, lat, stat_code))

    #or use the convenience function
    #write_statsll(loc, fname, event_stats1)
    write_stats_kml(loc, fname_statsll.split(".")[0]+".kml", event_stats1)

    std_out.close()
    final_time=time()
    print("Done in {:.1f} secs\n".format(final_time - init_time))

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

    #Only LOC needs to be changed in this example, LOC is where the geoNet data
    #files are placed after download
    #LOC="/".join([os.getcwd(), "Vol1", "data"])

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
            #record might then be empty, skipp processing in such a case
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
        stat_code = station_file_name.split(".")[0].split("_")[2]

        try:
            #expects that velBB etc directories already exist, created when getData.py is used
            pgf.save2disk(LOC+"/velBB/", stat_code, 'velBB')
            pgf.save2disk(LOC+"/velLF/", stat_code, 'velLF')
            pgf.save2disk(LOC+"/accBB/", stat_code, 'accBB')

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

    #import numpy as np
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
    from geoNet.utils import get_sorted_stats_code,  read_statsll, get_processed_stats_list

    init_time = time()
    #parent_dir_loc="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017/Vol1/data"
    #parent_dir_loc="/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105/Vol1/data"
    #plot_dir_accBB="accBB"
    #plot_dir_velBB="velBB"
    #plot_dir_accLF_0pt25="accLF_0pt25"
    #plot_dir_velLF_0pt25="velLF_0pt25"
    #plot_dir="velBB"
    loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])
    #loc_accLF_0pt25="/".join([parent_dir_loc, plot_dir_accLF_0pt25])
    #loc_velLF_0pt25="/".join([parent_dir_loc, plot_dir_velLF_0pt25])

    loc_acc = loc_accBB
    loc_vel = loc_velBB


    #stats_dict = read_statsll("/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105",
    #                          'geonet_stations_20170117.ll')
    stats_dict = read_statsll(loc_statsll, fname_statsll)
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)

    #stations are sorted according to PGV
    sorted_stats_code = get_sorted_stats_code(loc_velBB,stats_dict,comp='geom')
    pdf_acc = PdfPages('plots_acc.pdf')
    pdf_vel = PdfPages('plots_vel.pdf')
    for stat in sorted_stats_code:
        #stat = stat.strip("\n")
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

    #parent_dir_loc="/nesi/projects/nesi00213/ObservedGroundMotions/Mw5pt95_20150105_174841/Vol1/data"
    #parent_dir_loc="/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105/Vol1/data"
    ##plot_dir="velLF"
    #plot_dir_accBB="accBB"
    #plot_dir_velBB="velBB"
    #plot_dir_accLF_0pt25="accLF_0pt25"
    #plot_dir="velBB"
    loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])
    #loc_accLF_0pt25="/".join([parent_dir_loc, plot_dir_accLF_0pt25])

    loc_acc = loc_accBB


    #stats_dict = read_statsll("/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105",
    #                          'geonet_stations_20170117.ll')
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
    with open("event_obs_PGVs.txt", 'w') as f:
        for stat_data in event_velBB:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            ext_vel_000, _ = get_extremum(stat_data['000'])
            ext_vel_090, _ = get_extremum(stat_data['090'])
            max_vel_geom = np.sqrt(np.abs(ext_vel_000*ext_vel_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line+= " {:^15.6f}".format(max_vel_geom)
            f.write(line+"\n")

    #write PGAs
    event_accBB = get_event_data(loc_accBB, sorted_stats_code)
    with open("event_obs_PGAs.txt", 'w') as f:
        for stat_data in event_accBB:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            ext_acc_000, _ = get_extremum(stat_data['000'])
            ext_acc_090, _ = get_extremum(stat_data['090'])
            max_acc_geom = np.sqrt(np.abs(ext_acc_000*ext_acc_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line+= " {:^15.6f}".format(max_acc_geom)
            f.write(line+"\n")

    #sys.exit("getting out")


    #need new iterator object
    #write PSAs
    period=np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])
    event_PSA   = get_event_PSA(get_event_data(loc_accBB, sorted_stats_code),
                                period, xi=0.05, m=1., gamma=0.5, beta=0.25)


    with open("event_obs_PSAs.txt", 'w') as f:
        f.write((len(period)*" {:5.1f} ").format(*period))
        f.write("\n")
        for stat_data in event_PSA:
            stat_code=stat_data['name']
            lon, lat = stats_dict[stat_code]
            pSA_geom = stat_data['geom']
            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line+= (len(period)*" {:^15.6f} ").format(*pSA_geom)
            f.write(line+"\n")

    final_time = time()
    print("Done in {:.1f} secs.\n".format(final_time - init_time))
    return

if __name__ == "__main__":

    loc = os.getcwd()
    #BASE_URL="ftp://ftp.geonet.org.nz/strong/processed/Proc/2017/02_Feb/2017-02-01_102129/Vol1/data/"
    BASE_URL="ftp://ftp.geonet.org.nz/strong/processed/Proc/2017/02_Feb/2017-02-02_074142/Vol1/data/"
    getData(loc, BASE_URL)

    fname_statsll="event_stats.ll"
    loc_statsll=os.getcwd()
    #loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo"
    loc_all_geoNet_stats="."
    fname_all_geoNet_stats="all_geoNet_stats_2016-12-20.ll"
    loc_V1A="/".join([os.getcwd(), "Vol1", "data"])

    event_statsll(fname_statsll, loc_statsll,
                  loc_all_geoNet_stats, fname_all_geoNet_stats,
                  loc_V1A)


    #LOC="/".join([os.getcwd(), "Vol1", "data"])
    processData(loc_V1A)

    plot_accvel(loc_V1A, "accBB", "velBB",
                loc_statsll, fname_statsll)

    plot_psa(loc_V1A, "accBB", "velBB",
             loc_statsll, fname_statsll)

    IMsOnMap(loc_V1A, "accBB", "velBB",
             loc_statsll, fname_statsll)



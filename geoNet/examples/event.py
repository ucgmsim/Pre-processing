"""
Use in analysis of event to produce data files for plotting etc
"""

import numpy as np
import os
from time import time
from glob import glob
import scipy
import sys
import operator as op

# geoNet imports
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
    init_time = time()
    print("Downloading data ...")

    sg.make_dataPlot_dirs_v2(loc=loc, num_vol=1)

    std_out, std_err = sg.get_geoNet_data(
        loc="/".join([loc, "Vol1"]),
        geoNet_dir="data",
        geoNet_url=BASE_URL,
        wget_options="-N",
    )

    with open("getData_stdout.txt", "w") as f:
        f.writelines(std_out)

    with open("getData_stderr.txt", "w") as f:
        f.writelines(std_err)

    final_time = time()
    print("Finished downloading data in {:.1f} secs\n".format(final_time - init_time))

    return


def event_statsll(fname, loc, loc_all_geoNet_stats, fname_all_geoNet_stats, loc_V1A):
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
    init_time = time()
    print("Creating stations list ...")

    std_out = open("event_statsll_stdout.txt", "w")
    # Firs get statsll dictionary by assigning lon, lat from the saved list known
    # to be in WGS84 coordinates
    (event_stats1, fname_statsll) = sg.statsll_from_V1A_v2(
        loc_all_geoNet_stats,
        fname_all_geoNet_stats,
        loc_V1A,
        save_stats=False,
        fname=None,
        loc=loc,
    )

    # Now create another dictionary in which the lon, lat are read from .V1A files
    # which may or may not be in WGS84 coordinates. Unless something goes wrong
    # this should contain all the downloaded stations
    # Note fname=None crashed
    event_stats2, _ = sg.statsll_from_V1A(
        loc_V1A, save_stats=False, fname=None, loc=loc
    )

    # Find stats that are in event_stat2 but not in event_stats1
    stat_codes1 = set(event_stats1.keys())
    stat_codes2 = set(event_stats2.keys())

    # perform check
    if not stat_codes1.issubset(stat_codes2):
        std_out.write("Some station (lon, lat) were not read from .V1A files\n")

    for stat_code in stat_codes2 - stat_codes1:
        event_stats1[stat_code] = event_stats2[stat_code]
    # write statsll for this event
    with open("/".join([loc, fname]), "w") as f:
        for stat_code in event_stats1:
            (lon, lat) = event_stats1[stat_code]
            f.write("{:<15.4f} {:^15.4f} {:^10s} \n".format(lon, lat, stat_code))

    # write kml file to view with google earth/map
    write_stats_kml(loc, fname_statsll.split(".")[0] + ".kml", event_stats1)

    std_out.close()
    final_time = time()
    print("Done in {:.1f} secs\n".format(final_time - init_time))

    return


def processData(LOC):
    """
    LOC:
        where .V1A data files are saved in getData.
    """
    init_time = time()
    print("Processing .V1A files ...")

    std_err = open("processData_stderr.txt", "w")
    std_out = open("processData_stdout.txt", "w")
    sys.stdout = std_out

    FILE_NAMES = []
    event_stats_V1A = glob("/".join([LOC, "*.V1A"]))
    FILE_NAMES = [os.path.basename(_) for _ in event_stats_V1A]

    std_out.write("\n Processing %d stations in Vol1 data ...\n" % len(FILE_NAMES))

    for station_file_name in FILE_NAMES:
        std_out.write("\n**************************\n")
        std_out.write("%s" % station_file_name)
        std_out.write("\n**************************\n")
        try:
            gf = GeoNet_File(station_file_name, LOC, vol=1)
            if gf.comp_1st.acc.size < 5.0 / gf.comp_1st.delta_t:
                std_out.write("%s has less than 5 secs of data\n" % station_file_name)
                std_out.write("skipping %s\n" % station_file_name)
                continue

            # When appended zeroes at the beginning of the record are removed, the
            # record might then be empty, skipp processing in such a case.
            # The number 10 is rather arbirary, 1 should suffice
            agf = adjust_gf_for_time_delay(gf)
            if agf.comp_1st.acc.size <= 1:
                std_out.write("no elements in %s. Skipping it.\n" % station_file_name)
                continue

            gf.comp_1st.acc = scipy.signal.detrend(gf.comp_1st.acc, type="linear")
            gf.comp_2nd.acc = scipy.signal.detrend(gf.comp_2nd.acc, type="linear")
            gf.comp_up.acc = scipy.signal.detrend(gf.comp_up.acc, type="linear")

            gf.comp_1st.acc -= gf.comp_1st.acc.mean()
            gf.comp_2nd.acc -= gf.comp_2nd.acc.mean()
            gf.comp_up.acc -= gf.comp_up.acc.mean()

            pgf = Process(gf, lowcut=0.05, gpInt=False)
            # pgf = Process(gf, lowcut=0.05, ft=0.25)
        except Exception as e:
            std_err.write(str(e))
            std_err.write("%s is problematic, skipping it.\n" % station_file_name)
            continue
            # raise

        # extract stat code
        # stat_code = station_file_name.split(".")[0].split("_")[2]
        # This could be placed in a try, except block
        YYYYMMDD, HHMM, stat_code = utils.parse_V1A_fname(station_file_name)

        try:
            # expects that velBB etc directories already exist, created when getData.py is used
            pgf.save2disk(LOC + "/velBB/", stat_code, "velBB")
            pgf.save2disk(LOC + "/velLF/", stat_code, "velLF")
            pgf.save2disk(LOC + "/accBB/", stat_code, "accBB")

        except Exception as e:
            std_err.write(str(e))
            std_err.write("Skipping this station %s\n" % station_file_name)
            continue

    std_out.close()
    std_err.close()

    sys.stdout = sys.__stdout__
    final_time = time()
    print("Done processing in {:.1f} secs\n".format(final_time - init_time))

    return


def plot_accvel(
    parent_dir_loc, plot_dir_accBB, plot_dir_velBB, loc_statsll, fname_statsll
):
    """
    Plots all the SMS data in a given directory by reading all *.000 files and
    sorting the SMSs in order of increasing PGV
    """
    init_time = time()
    print("Creating acc and vel plots for all Observed SMSs ...")

    import matplotlib

    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use("Agg")
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
    from geoNet.utils import (
        get_sorted_stats_code,
        read_statsll,
        get_processed_stats_list,
    )

    init_time = time()
    loc_accBB = "/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB = "/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB
    loc_vel = loc_velBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    # Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)

    # stations are sorted according to PGV
    sorted_stats_code = get_sorted_stats_code(loc_velBB, stats_dict, comp="geom")
    pdf_acc = PdfPages("plots_acc.pdf")
    pdf_vel = PdfPages("plots_vel.pdf")
    for stat in sorted_stats_code:
        stat_acc, fig_acc, ax_acc = get_stat_acc_plot(loc_acc, stat["name"])
        pdf_acc.savefig()

        stat_acc, fig_vel, ax_vel = get_stat_vel_plot(loc_vel, stat["name"])
        pdf_vel.savefig()
        plt.close("all")

    pdf_acc.close()
    pdf_vel.close()

    final_time = time()
    print("Done plotting acc, vel in {:.1f} secs.\n".format(final_time - init_time))
    return


def plot_psa(
    parent_dir_loc, plot_dir_accBB, plot_dir_velBB, loc_statsll, fname_statsll
):
    """ """
    import matplotlib

    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use("Agg")
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from geoNet.putils import get_stat_PSA_plot
    from geoNet.utils import (
        get_sorted_stats_code,
        read_statsll,
        get_processed_stats_list,
    )

    init_time = time()
    print("Creating pSA plots for all Observed SMSs ...")

    loc_accBB = "/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB = "/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    # Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)
    # stations are sorted according to PGV
    sorted_stats_code = get_sorted_stats_code(loc_velBB, stats_dict)
    pdf_psa = PdfPages("plots_psa.pdf")
    for stat in sorted_stats_code:
        stat, PSA, fig_psa, ax = get_stat_PSA_plot(loc_acc, stat["name"])
        pdf_psa.savefig()

        plt.close("all")

    pdf_psa.close()

    final_time = time()
    print("Done pSA plots in {:.1f} secs.\n".format(final_time - init_time))

    return


def IMsOnMap(
    parent_dir_loc, plot_dir_accBB, plot_dir_velBB, loc_statsll, fname_statsll
):
    """
    Writes PGV, PGA, pSA intensity measures to files used in GMT plotting scripts
    """

    from geoNet.utils import (
        get_sorted_stats_code,
        read_statsll,
        get_event_data,
        get_event_PSA,
        get_extremum,
        get_processed_stats_list,
    )

    init_time = time()
    print("Writting IMs for GMT plotting ...")

    loc_accBB = "/".join([parent_dir_loc, plot_dir_accBB])
    loc_velBB = "/".join([parent_dir_loc, plot_dir_velBB])

    loc_acc = loc_accBB
    loc_vel = loc_velBB

    stats_dict = read_statsll(loc_statsll, fname_statsll)
    # Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict = get_processed_stats_list(loc_velBB, stats_dict, verbose=True)

    # stations are sorted according to PGV
    # Altough not required helps identify maximum and minimum IM interval
    sorted_stats_code = get_sorted_stats_code(loc_velBB, stats_dict, comp="geom")
    sorted_stats_code = [_["name"] for _ in sorted_stats_code]

    # write PGVs
    event_velBB = get_event_data(loc_velBB, sorted_stats_code)
    with open("event_obs_PGVs.txt", "w") as f:
        for stat_data in event_velBB:
            stat_code = stat_data["name"]
            lon, lat = stats_dict[stat_code]
            ext_vel_000, _ = get_extremum(stat_data["000"])
            ext_vel_090, _ = get_extremum(stat_data["090"])
            max_vel_geom = np.sqrt(np.abs(ext_vel_000 * ext_vel_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line += " {:^15.6f}".format(max_vel_geom)
            f.write(line + "\n")

    # write PGAs
    event_accBB = get_event_data(loc_accBB, sorted_stats_code)
    with open("event_obs_PGAs.txt", "w") as f:
        for stat_data in event_accBB:
            stat_code = stat_data["name"]
            lon, lat = stats_dict[stat_code]
            ext_acc_000, _ = get_extremum(stat_data["000"])
            ext_acc_090, _ = get_extremum(stat_data["090"])
            max_acc_geom = np.sqrt(np.abs(ext_acc_000 * ext_acc_090))

            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line += " {:^15.6f}".format(max_acc_geom)
            f.write(line + "\n")

    # sys.exit("getting out")

    # need new iterator object
    # write PSAs
    period = np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])
    event_PSA = get_event_PSA(
        get_event_data(loc_accBB, sorted_stats_code),
        period,
        xi=0.05,
        m=1.0,
        gamma=0.5,
        beta=0.25,
    )

    with open("event_obs_PSAs.txt", "w") as f:
        f.write((len(period) * " {:5.1f} ").format(*period))
        f.write("\n")
        for stat_data in event_PSA:
            stat_code = stat_data["name"]
            lon, lat = stats_dict[stat_code]
            pSA_geom = stat_data["geom"]
            line = "{:10.4f} {:10.4f}".format(lon, lat)
            line += (len(period) * " {:^15.6f} ").format(*pSA_geom)
            f.write(line + "\n")

    final_time = time()
    print("Done in {:.1f} secs.\n".format(final_time - init_time))
    return


def IMsOnMap_ratios(
    parent_dir_loc_obs,
    plot_dir_accBB_obs,
    plot_dir_velBB_obs,
    loc_statsll_obs,
    fname_statsll_obs,
    parent_dir_loc_sim,
    plot_dir_accBB_sim,
    plot_dir_velBB_sim,
    loc_statsll_sim,
    fname_statsll_sim,
):
    """
    Writes PGV, PGA, pSA intensity measures to files used in GMT plotting scripts
    """

    from geoNet.utils import (
        get_sorted_stats_code,
        read_statsll,
        get_event_data,
        get_event_PSA,
        get_extremum,
        get_processed_stats_list,
        get_stat_data,
        get_PSA,
    )

    init_time = time()
    print("Writting IMs for GMT plotting ...")

    loc_accBB_obs = "/".join([parent_dir_loc_obs, plot_dir_accBB_obs])
    loc_velBB_obs = "/".join([parent_dir_loc_obs, plot_dir_velBB_obs])

    loc_acc_obs = loc_accBB_obs
    loc_vel_obs = loc_velBB_obs

    loc_accBB_sim = "/".join([parent_dir_loc_sim, plot_dir_accBB_sim])
    loc_velBB_sim = "/".join([parent_dir_loc_sim, plot_dir_velBB_sim])

    loc_acc_sim = loc_accBB_sim
    loc_vel_sim = loc_velBB_sim

    stats_dict_sim = read_statsll(loc_statsll_sim, fname_statsll_sim)
    stats_dict_obs = read_statsll(loc_statsll_obs, fname_statsll_obs)
    # Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict_obs = get_processed_stats_list(
        loc_velBB_obs, stats_dict_obs, verbose=True
    )

    # stations are sorted according to PGV
    # Altough not required helps identify maximum and minimum IM interval
    sorted_stats_code = get_sorted_stats_code(
        loc_velBB_obs, stats_dict_obs, comp="geom"
    )
    sorted_stats_code = [_["name"] for _ in sorted_stats_code]

    # write PGVs
    event_velBB = get_event_data(loc_velBB_obs, sorted_stats_code)
    with open("event_ratios_PGVs.txt", "w") as f:
        lines = []
        for stat_data in event_velBB:
            stat_code = stat_data["name"]
            # Only take ratios if SMS is both in sims and obs
            if not stats_dict_sim.has_key(stat_code):
                continue
            lon, lat = stats_dict_obs[stat_code]
            ext_vel_000, _ = get_extremum(stat_data["000"])
            ext_vel_090, _ = get_extremum(stat_data["090"])
            max_vel_geom = np.sqrt(np.abs(ext_vel_000 * ext_vel_090))

            # simulations
            stat_data_sim = get_stat_data(loc_vel_sim, stat_code)
            ext_vel_000_sim, _ = get_extremum(stat_data_sim["000"])
            ext_vel_090_sim, _ = get_extremum(stat_data_sim["090"])
            max_vel_geom_sim = np.sqrt(np.abs(ext_vel_000_sim * ext_vel_090_sim))
            lines.append([lon, lat, np.log(max_vel_geom / max_vel_geom_sim)])

        lines = sorted(lines, key=op.itemgetter(2))
        for line in lines:
            # line = "{:10.4f} {:10.4f}".format(lon, lat)
            ##line+= " {:^15.6f}".format(max_vel_geom)
            # line+= " {:^15.6f}".format(np.log(max_vel_geom/max_vel_geom_sim))
            line = "{:10.4f} {:10.4f} {:^15.6f}".format(*line)
            f.write(line + "\n")

    # write PGAs
    event_accBB = get_event_data(loc_accBB_obs, sorted_stats_code)
    with open("event_ratios_PGAs.txt", "w") as f:
        lines = []
        for stat_data in event_accBB:
            stat_code = stat_data["name"]
            # Only take ratios if SMS is both in sims and obs
            if stats_dict_sim.has_key(stat_code) is False:
                continue
            lon, lat = stats_dict_obs[stat_code]
            ext_acc_000, _ = get_extremum(stat_data["000"])
            ext_acc_090, _ = get_extremum(stat_data["090"])
            max_acc_geom = np.sqrt(np.abs(ext_acc_000 * ext_acc_090))

            # simulations
            stat_data_sim = get_stat_data(loc_acc_sim, stat_code)
            ext_acc_000_sim, _ = get_extremum(stat_data_sim["000"])
            ext_acc_090_sim, _ = get_extremum(stat_data_sim["090"])
            max_acc_geom_sim = np.sqrt(np.abs(ext_acc_000_sim * ext_acc_090_sim))
            # simulation acc are in cm/s^2 change to g units
            max_acc_geom_sim /= 981.0

            lines.append([lon, lat, np.log(max_acc_geom / max_acc_geom_sim)])

        lines = sorted(lines, key=op.itemgetter(2))
        for line in lines:
            # line = "{:10.4f} {:10.4f}".format(lon, lat)
            # line+= " {:^15.6f}".format(np.log(max_acc_geom/max_acc_geom_sim))
            line = "{:10.4f} {:10.4f} {:^15.6f}".format(*line)
            f.write(line + "\n")

    # need new iterator object
    # write PSAs
    period = np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])
    event_PSA = get_event_PSA(
        get_event_data(loc_accBB_obs, sorted_stats_code),
        period,
        xi=0.05,
        m=1.0,
        gamma=0.5,
        beta=0.25,
    )

    with open("event_ratios_PSAs.txt", "w") as f:
        lines = []
        f.write((len(period) * " {:5.1f} ").format(*period))
        f.write("\n")
        for stat_data in event_PSA:
            stat_code = stat_data["name"]
            # Only take ratios if SMS is both in sims and obs
            if stats_dict_sim.has_key(stat_code) is False:
                continue
            lon, lat = stats_dict_obs[stat_code]
            pSA_geom = stat_data["geom"]

            # simulations
            stat_data_sim = get_stat_data(loc_acc_sim, stat_code)
            # simulation acc are in cm/s^2 change to g units
            for key in ["000", "090", "ver"]:
                stat_data_sim[key] /= 981.0
            dt = stat_data_sim["t"][1] - stat_data_sim["t"][0]
            stat_PSA_sim = get_PSA(stat_data_sim, dt, period)
            pSA_geom_sim = stat_PSA_sim["geom"]

            lines.append([lon, lat] + list(np.log(pSA_geom / pSA_geom_sim)))

        # sorted pSA according to smallest period
        lines = sorted(lines, key=op.itemgetter(2))
        for line in lines:
            # line = "{:10.4f} {:10.4f}".format(lon, lat)
            # line+= (len(period)*" {:^15.6f} ").format(*(np.log(pSA_geom/pSA_geom_sim)))
            line = ("{:10.4f} {:10.4f}" + len(period) * " {:^15.6f}").format(*line)
            f.write(line + "\n")

    final_time = time()
    print("Done in {:.1f} secs.\n".format(final_time - init_time))
    return


def keyValueFromTxt(fname):
    """
    Parses file that has the form key=value and returns a dictionary
    """
    keyValue = dict()
    fname = os.path.abspath(fname)
    print("Reading input from {:s}\n".format(fname))
    with open(fname, "r") as f:
        for line in f:
            # remove white spaces
            if line.startswith("#") or line.startswith("%"):
                continue
            if line in ["\n", "\r", "\rn"]:
                continue
            line = line.strip()
            line = line.replace(" ", "")
            line = line.replace('"', "")
            key, value = line.split("=")
            keyValue[key] = value

    return keyValue


def run(arg):
    """
    main function that runs realtime.py
    """
    keyValue = keyValueFromTxt(arg)

    loc = keyValue["loc"]
    BASE_URL = keyValue["BASE_URL"]
    loc_statsll = keyValue["loc_statsll"]
    fname_statsll = keyValue["fname_statsll"]
    loc_all_geoNet_stats = keyValue["loc_all_geoNet_stats"]
    fname_all_geoNet_stats = keyValue["fname_all_geoNet_stats"]
    loc_V1A = keyValue["loc_V1A"]
    obs_velDir = keyValue["obs_velDir"]
    obs_accDir = keyValue["obs_accDir"]

    sim_velDir = keyValue["sim_velDir"]
    sim_accDir = keyValue["sim_accDir"]
    lf_sim_dir = keyValue["lf_sim_dir"]
    hf_dir = keyValue["hf_dir"]
    bb_dir = keyValue["bb_dir"]
    loc_statsll_sim = keyValue["loc_statsll_sim"]
    fname_statsll_sim = keyValue["fname_statsll_sim"]

    srf_dir = keyValue["srf_dir"]
    srf_file = keyValue["srf_file"]
    # event_statsll(fname_statsll, loc_statsll,
    #              loc_all_geoNet_stats, fname_all_geoNet_stats,
    #              loc_V1A)

    # processData(loc_V1A)

    # plot_accvel(loc_V1A, obs_accDir, obs_velDir,
    #            loc_statsll, fname_statsll)

    # plot_psa(loc_V1A, obs_accDir, obs_velDir,
    #         loc_statsll, fname_statsll)

    # IMsOnMap(loc_V1A, obs_accDir, obs_velDir,
    #         loc_statsll, fname_statsll)

    IMsOnMap_ratios(
        loc_V1A,
        obs_accDir,
        obs_velDir,
        loc_statsll,
        fname_statsll,
        bb_dir,
        sim_accDir,
        sim_velDir,
        loc_statsll_sim,
        fname_statsll_sim,
    )

    return


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Use for event analyses.")
    parser.add_argument(
        "-f", "--fname", help="name of input file is required", required=True
    )
    args = vars(parser.parse_args())
    run(args["fname"])

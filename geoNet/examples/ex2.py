"""
Calculate and plot Bias= log(pSA_obs/pSA_sim).
Similar to ex1.py but sorts stats according to Rrup
"""
import sys
import os
import numpy as np
import glob
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")
from matplotlib import pylab as plt

from geoNet.utils import (
    read_statsll,
    get_processed_stats_list,
    get_SMS_Rrups,
    get_SMS_pSA,
    get_bias,
)
from geoNet.putils import plot_bias
from geoNet.gmpe import readStationFile as rsf
from geoNet.gmpe.calculateGMPE import set_faultprop


# *****************************************************************************
#                       OBSERVATIONS
# *****************************************************************************
parent_dir_loc_obs = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt6_20100904_103801_11Jan2017/Vol1/data"

stats_dict_obs = read_statsll(
    "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017",
    "20100904_103801_eventStats_2017-01-10.ll",
)


plot_dir_accBB_obs = "accBB_1"
plot_dir_velBB_obs = "velBB_1"

loc_accBB_obs = "/".join([parent_dir_loc_obs, plot_dir_accBB_obs])
loc_velBB_obs = "/".join([parent_dir_loc_obs, plot_dir_velBB_obs])

loc_acc_obs = loc_accBB_obs
loc_vel_obs = loc_velBB_obs

# reset stats_dict to those processed
stats_dict_obs = get_processed_stats_list(loc_vel_obs, stats_dict_obs)


# *****************************************************************************
#                       SIMULATIONS
# *****************************************************************************

base_dir_sim = "/nesi/projects/nesi00213/RunFolder/ahsan"
parent_dir_loc_sim = "/".join(
    [base_dir_sim, "ptSource_2010Sep04_v2_m4pt6_VMv1p64_FVM-h0p100_170123"]
)

stats_dict_sim = read_statsll(parent_dir_loc_sim, "HF_BB_stats.ll")


plot_dir_accBB_sim = "Acc"
plot_dir_velBB_sim = "Vel"

loc_accBB_sim = glob.glob(
    "/".join([parent_dir_loc_sim, "BB", "*", plot_dir_accBB_sim])
)[0]
loc_velBB_sim = glob.glob(
    "/".join([parent_dir_loc_sim, "BB", "*", plot_dir_velBB_sim])
)[0]

loc_acc_sim = loc_accBB_sim
loc_vel_sim = loc_velBB_sim


# *****************************************************************************
#                   Perform calculations and plot
#
# *****************************************************************************
# (1) construct a dictionary that only contains SMS present both in simulations and observations
stats_codes_simObs = list(set(stats_dict_obs) & set(stats_dict_sim))
stats_dict_simObs = {}
for stat_code in stats_codes_simObs:
    lon, lat = stats_dict_sim[stat_code]
    stats_dict_simObs[stat_code] = (lon, lat)

# (2) Sort the SMSs according to distance, need to know fault/source properties
srf_fname = "/nesi/projects/nesi00213/RupModel/ahsan/2010Sep04_m4pt6/ptSource/Srf/2010Sep04_v2_m4pt6.srf"
FiniteFault = rsf.readSrfFile(srf_fname)
FiniteFault = rsf.Points_np(FiniteFault)
# set fault properties
faultprop = set_faultprop(Mw=4.6, rake=173.0, dip=74.0, Ztor=8.0)

sorted_stat_codes, sms_rrups = get_SMS_Rrups(stats_dict_simObs, FiniteFault)

# (3) calculate pSA for simulations and observation for the above SMSs
# Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
# one or the other
g = 981.0  # cm/s^2
periods = np.logspace(start=np.log10(0.01), stop=np.log10(10.0), num=100, base=10)

pSA_sim = get_SMS_pSA(sorted_stat_codes, periods, loc_acc_sim, comp="geom") / g
pSA_obs = get_SMS_pSA(sorted_stat_codes, periods, loc_acc_obs, comp="geom")
bias, std = get_bias(pSA_obs, pSA_sim, rescale=False)

# (4) Finally plot the bias
figBias, axBias = plot_bias(
    periods, bias, std, savefig=True, figName="ex2_bias", ext="png"
)

import sys
import os
import numpy as np
import glob
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

from geoNet.utils import read_statsll, get_processed_stats_list
from geoNet.gmpe import readStationFile as rsf
from geoNet.gmpe.calculateGMPE import set_faultprop 


#*****************************************************************************
#                       OBSERVATIONS
#*****************************************************************************
parent_dir_loc_obs="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt6_20100904_103801_11Jan2017/Vol1/data"

stats_dict_obs = read_statsll("/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017",
"20100904_103801_eventStats_2017-01-10.ll")


plot_dir_accBB_obs="accBB_1"
plot_dir_velBB_obs="velBB_1"

loc_accBB_obs="/".join([parent_dir_loc_obs, plot_dir_accBB_obs])
loc_velBB_obs="/".join([parent_dir_loc_obs, plot_dir_velBB_obs])

loc_acc_obs = loc_accBB_obs
loc_vel_obs = loc_velBB_obs

#reset stats_dict to those processed
stats_dict_obs=get_processed_stats_list(loc_vel_obs,stats_dict_obs)



#*****************************************************************************
#                       SIMULATIONS
#*****************************************************************************

base_dir_sim="/nesi/projects/nesi00213/RunFolder/ahsan"
parent_dir_loc_sim="/".join([base_dir_sim, "ptSource_2010Sep04_v2_m4pt6_VMv1p64_FVM-h0p100_170123"])

stats_dict_sim = read_statsll(parent_dir_loc_sim, "HF_BB_stats.ll")


plot_dir_accBB_sim="Acc"
plot_dir_velBB_sim="Vel"

loc_accBB_sim=glob.glob("/".join([parent_dir_loc_sim, "BB", "*", plot_dir_accBB_sim]))[0]
loc_velBB_sim=glob.glob("/".join([parent_dir_loc_sim, "BB", "*", plot_dir_velBB_sim]))[0]

loc_acc_sim = loc_accBB_sim
loc_vel_sim = loc_velBB_sim

#*****************************************************************************

#*****************************************************************************
#                   Read the standard rupture file
#*****************************************************************************
srf_fname="/nesi/projects/nesi00213/RupModel/ahsan/2010Sep04_m4pt6/ptSource/Srf/2010Sep04_v2_m4pt6.srf"
FiniteFault = rsf.readSrfFile(srf_fname)
FiniteFault = rsf.Points_np(FiniteFault)

#set fault properties
faultprop = set_faultprop(Mw=4.6, rake=173., dip=74., Ztor=8.)


periods_gmpe=np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])
Rrups_gmpe = np.logspace(np.log10(5.),np.log10(100.),30)

#get_empIM_v2(Rrup, period, faultprop, Rjb=None, Rtvz=0., V30measured=0., V30=250.)
empIM_values, empIM_sigmas =  get_empIM_v2(Rrups_gmpe, periods_gmpe, faultprop)

#def plot_IMvsRrup(Rrup, IM)

fig, ax = plt.subplots()
ax.plot(Rrups_gmpe, empIM_values[:,0],c='k', ls='solid')
ax.plot(Rrups_gmpe, empIM_values[:,0]*np.exp(empIM_sigmas[:,0]), c='k',ls='dashed')
ax.plot(Rrups_gmpe, empIM_values[:,0]*np.exp(-empIM_sigmas[:,0]), c='k', ls='dashed')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(Rrups_gmpe[0],Rrups_gmpe[-1])
fig.savefig("gmpe_pSA0.1.png")

stats_codes_simObs = list(set(stats_dict) & set(stats_dict_sim))
stats_dict_simObs = {}
for stat_code in stats_codes_simObs:
    lon, lat = stats_dict_sim[stat_code]
    stats_dict_simObs[stat_code] = (lon, lat)

sorted_stat_codes, sms_rrups = get_SMS_Rrups(stats_dict_simObs, FiniteFault)
period=periods_gmpe[0]
g=981. #cm/s^2
pSA_sim = get_SMS_pSA(sorted_stat_codes, period, loc_acc_sim, comp='geom')/g
pSA_obs = get_SMS_pSA(sorted_stat_codes, period, loc_acc, comp='geom')
fig, ax = plot_SMS_IM(sms_rrups[:,0], pSA_sim[:,0], pSA_obs[:,0])

#get_empIM_v2(Rrup, period, faultprop, Rjb=None, Rtvz=0., V30measured=0., V30=250.)
#empIM_values, empIM_sigmas =  get_empIM_v2(Rrups_gmpe, periods_gmpe, faultprop)
pSA_gmpe, pSA_gmpe_std = get_empIM_v2(Rrups_gmpe, period, faultprop)
fig, ax = plot_IMvsRrup(Rrups_gmpe, pSA_gmpe, pSA_gmpe_std, fig=fig, ax=ax)
ax.legend(loc="best", scatterpoints=1)
fig.savefig("pSA_ex.png")

fig, ax = plot_SMS_IM_ratio(sms_rrups[:,0], pSA_obs[:,0]/pSA_sim[:,0])
fig.savefig("pSA_ex_ratio.png")
import sys
sys.exit("outa here")
periods = np.logspace(start=np.log10(0.01), stop=np.log10(10.), num=100, base=10)
#periods = np.concatenate([
#          np.logspace(start=np.log10(0.01), stop=np.log10(0.1), num=100, base=10, endpoint=False),
#          np.logspace(start=np.log10(0.1), stop=np.log10(1.), num=100, base=10, endpoint=False),
#          np.logspace(start=np.log10(1.), stop=np.log10(10.), num=100, base=10, endpoint=True)
#          ])
#periods = np.concatenate([
#          np.logspace(start=np.log10(0.01), stop=np.log10(0.1), num=33, base=10, endpoint=False),
#          np.logspace(start=np.log10(0.1), stop=np.log10(1.), num=34, base=10, endpoint=False),
#          np.logspace(start=np.log10(1.), stop=np.log10(10.), num=33, base=10, endpoint=True)
#          ])
pSA_sim = get_SMS_pSA(sorted_stat_codes, periods, loc_acc_sim, comp='geom')
pSA_obs = get_SMS_pSA(sorted_stat_codes, periods, loc_acc, comp='geom')
bias, std = get_bias(pSA_obs, pSA_sim)
figBias, axBias = plot_bias(bias, std, savefig=True)

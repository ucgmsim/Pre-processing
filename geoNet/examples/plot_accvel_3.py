import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
from time import time

from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
from geoNet.utils import get_sorted_stats_code,  read_statsll

init_time = time()
parent_dir_loc="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw7pt5_20161113_110256/Vol1/data"
plot_dir_accBB="accBB"
plot_dir_velBB="velBB"
loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])

loc_acc = loc_accBB
loc_vel = loc_velBB

#for individual plots use example below
stat_POTS, fig_acc_POTS, ax_acc_POTS = get_stat_acc_plot(loc_acc, "POTS") #"WDFS", MTHS
ax_acc_POTS[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_POTS[0].set_xlim(50,200.)
ax_acc_POTS[1].set_xlim(50,200.)
ax_acc_POTS[2].set_xlim(50,200.)
##customize y-limits
ax_acc_POTS[0].set_ylim(-0.2,0.2)
ax_acc_POTS[1].set_ylim(-0.2,0.2)
ax_acc_POTS[2].set_ylim(-0.2,0.2)
fig_acc_POTS.savefig("acc_POTS.pdf")


stat_TEPS, fig_acc_TEPS, ax_acc_TEPS= get_stat_acc_plot(loc_acc, "TEPS") #"WDFS", MTHS
ax_acc_TEPS[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_TEPS[0].set_xlim(50.,200.)
ax_acc_TEPS[1].set_xlim(50,200.)
ax_acc_TEPS[2].set_xlim(50,200.)
##customize y-limits
ax_acc_TEPS[0].set_ylim(-0.20,0.20)
ax_acc_TEPS[1].set_ylim(-0.20,0.20)
ax_acc_TEPS[2].set_ylim(-0.2,0.2)
fig_acc_TEPS.savefig("acc_TEPS.pdf")


stat_WIGC, fig_acc_WIGC, ax_acc_WIGC = get_stat_acc_plot(loc_acc, "WIGC") #"WDFS", MTHS
ax_acc_WIGC[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_WIGC[0].set_xlim(0,150.)
ax_acc_WIGC[1].set_xlim(0,150.)
ax_acc_WIGC[2].set_xlim(0,150.)
#customize y-limits
ax_acc_WIGC[0].set_ylim(-1.1,1.1)
ax_acc_WIGC[1].set_ylim(-1.1,1.1)
ax_acc_WIGC[2].set_ylim(-1.1,1.1)
fig_acc_WIGC.savefig("acc_WIGC.pdf")


stat_CULC, fig_acc_CULC, ax_acc_CULC= get_stat_acc_plot(loc_acc, "CULC") #"WDFS", MTHS
ax_acc_CULC[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_CULC[0].set_xlim(0,150.)
ax_acc_CULC[1].set_xlim(0,150.)
ax_acc_CULC[2].set_xlim(0,150.)
#customize y-limits
ax_acc_CULC[0].set_ylim(-0.35,0.35)
ax_acc_CULC[1].set_ylim(-0.35,0.35)
ax_acc_CULC[2].set_ylim(-0.35,0.35)
fig_acc_CULC.savefig("acc_CULC.pdf")


stat_CECS, fig_acc_CECS, ax_acc_CECS= get_stat_acc_plot(loc_acc, "CECS") #"WDFS", MTHS
ax_acc_CECS[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_CECS[0].set_xlim(0,150.)
ax_acc_CECS[1].set_xlim(0,150.)
ax_acc_CECS[2].set_xlim(0,150.)
#customize y-limits
ax_acc_CECS[0].set_ylim(-0.35,0.35)
ax_acc_CECS[1].set_ylim(-0.35,0.35)
ax_acc_CECS[2].set_ylim(-0.15,0.15)
fig_acc_CECS.savefig("acc_CECS.pdf")


stat_HSES, fig_acc_HSES, ax_acc_HSES= get_stat_acc_plot(loc_acc, "HSES") #"WDFS", MTHS
ax_acc_HSES[2].set_xlabel('Time, t[s]')
#customize x-limits
ax_acc_HSES[0].set_xlim(0,150.)
ax_acc_HSES[1].set_xlim(0,150.)
ax_acc_HSES[2].set_xlim(0,150.)
#customize y-limits
ax_acc_HSES[0].set_ylim(-0.35,0.35)
ax_acc_HSES[1].set_ylim(-0.35,0.35)
ax_acc_HSES[2].set_ylim(-0.25,0.25)
fig_acc_HSES.savefig("acc_HSES.pdf")



#WTMC, KIKS, KEKS, WDFS
#stat_WTMC, fig_acc_WTMC, ax_acc_WTMC = get_stat_acc_plot(loc_acc, "WTMC") #"WDFS", MTHS
#ax_acc_WTMC[2].set_xlabel('Time, t[s]')
##customize x-limits
#ax_acc_WTMC[0].set_xlim(0,150.)
#ax_acc_WTMC[1].set_xlim(0,150.)
#ax_acc_WTMC[2].set_xlim(0,150.)
##customize y-limits
#ax_acc_WTMC[0].set_ylim(-1.5,1.5)
#ax_acc_WTMC[1].set_ylim(-1.5,1.5)
#ax_acc_WTMC[2].set_ylim(-3.5,3.5)
#fig_acc_WTMC.savefig("acc_WTMC.pdf")
#
#
#stat_KIKS, fig_acc_KIKS, ax_acc_KIKS = get_stat_acc_plot(loc_acc, "KIKS") #"WDFS", MTHS
#ax_acc_KIKS[2].set_xlabel('Time, t[s]')
##customize x-limits
#ax_acc_KIKS[0].set_xlim(0,150.)
#ax_acc_KIKS[1].set_xlim(0,150.)
#ax_acc_KIKS[2].set_xlim(0,150.)
##customize y-limits
#ax_acc_KIKS[0].set_ylim(-0.3,0.3)
#ax_acc_KIKS[1].set_ylim(-0.3,0.3)
#ax_acc_KIKS[2].set_ylim(-0.4,0.4)
#fig_acc_KIKS.savefig("acc_KIKS.pdf")
#
#
#stat_KEKS, fig_acc_KEKS, ax_acc_KEKS = get_stat_acc_plot(loc_acc, "KEKS") #"WDFS", MTHS
#ax_acc_KEKS[2].set_xlabel('Time, t[s]')
##customize x-limits
#ax_acc_KEKS[0].set_xlim(0,150.)
#ax_acc_KEKS[1].set_xlim(0,150.)
#ax_acc_KEKS[2].set_xlim(0,150.)
##customize y-limits
#ax_acc_KEKS[0].set_ylim(-1.6,1.6)
#ax_acc_KEKS[1].set_ylim(-1.6,1.6)
#ax_acc_KEKS[2].set_ylim(-0.5,0.5)
#fig_acc_KEKS.savefig("acc_KEKS.pdf")
#
#
#stat_WDFS, fig_acc_WDFS, ax_acc_WDFS = get_stat_acc_plot(loc_acc, "WDFS") #"WDFS", MTHS
#ax_acc_WDFS[2].set_xlabel('Time, t[s]')
##customize x-limits
#ax_acc_WDFS[0].set_xlim(0,150.)
#ax_acc_WDFS[1].set_xlim(0,150.)
#ax_acc_WDFS[2].set_xlim(0,150.)
##customize y-limits
#ax_acc_WDFS[0].set_ylim(-1.5,1.5)
#ax_acc_WDFS[1].set_ylim(-1.5,1.5)
#ax_acc_WDFS[2].set_ylim(-1.2,1.2)
#fig_acc_WDFS.savefig("acc_WDFS.pdf")
#plt.show()



#For velocity the call is to get_stat_vel_plot
#stat_WTMC, fig_vel_WTMC, ax_vel_WTMC = get_stat_vel_plot(loc_vel, "WTMC")
#fig_vel_WTMC.savefig("vel_MTHS.pdf")

plt.show()
final_time = time()
print("Done in {:.1f}".format(final_time - init_time))

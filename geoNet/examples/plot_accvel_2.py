import numpy as np
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
from time import time

from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
from geoNet.utils import get_sorted_stats_code, read_statsll

init_time = time()
parent_dir_loc = "/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw7pt5_20161113_110256/Vol1/data"
plot_dir_accBB = "accBB"
plot_dir_velBB = "velBB"
loc_accBB = "/".join([parent_dir_loc, plot_dir_accBB])
loc_velBB = "/".join([parent_dir_loc, plot_dir_velBB])

loc_acc = loc_accBB
loc_vel = loc_velBB

# for individual plots use example below
stat, fig_acc, ax = get_stat_acc_plot(loc_acc, "WTMC")  # "WDFS", MTHS
stat, fig_vel, ax = get_stat_vel_plot(loc_vel, "WTMC")

# fig_acc.savefig("acc_MTHS.pdf")
# fig_vel.savefig("vel_MTHS.pdf")

# final_time = time()
# print("Done in {:.1f}".format(final_time - init_time))

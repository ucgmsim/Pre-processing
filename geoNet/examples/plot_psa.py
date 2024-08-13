import numpy as np
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

from geoNet.putils import get_stat_PSA_plot
from geoNet.utils import get_sorted_stats_code, read_statsll

# parent_dir_loc="/nesi/projects/nesi00213/ObservedGroundMotions/Mw5pt95_20150105_174841/Vol1/data"
parent_dir_loc = (
    "/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105/Vol1/data"
)
# plot_dir="velLF"
plot_dir_accBB = "accBB"
plot_dir_velBB = "velBB"
# plot_dir_accLF_0pt25="accLF_0pt25"
# plot_dir="velBB"
loc_accBB = "/".join([parent_dir_loc, plot_dir_accBB])
loc_velBB = "/".join([parent_dir_loc, plot_dir_velBB])
# loc_accLF_0pt25="/".join([parent_dir_loc, plot_dir_accLF_0pt25])

loc_acc = loc_accBB

# for individual plots use example below
# stat, PSA, fig_psa, ax = get_stat_PSA_plot(loc_acc, "MTHS")
# fig_psa.savefig("psa")
# import sys
# sys.exit()
# To plot all stations in file all_stats use code below
# with open("all_stats.txt",'r') as f:
#    event_stats = f.readlines()

stats_dict = read_statsll(
    "/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105",
    "geonet_stations_20170117.ll",
)
# stations are sorted according to PSA
sorted_stats_code = get_sorted_stats_code(loc_velBB, stats_dict)
pdf_psa = PdfPages("plots_psa.pdf")
for stat in sorted_stats_code:
    stat, PSA, fig_psa, ax = get_stat_PSA_plot(loc_acc, stat["name"])
    pdf_psa.savefig()

    plt.close("all")

pdf_psa.close()

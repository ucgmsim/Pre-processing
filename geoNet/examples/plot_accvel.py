import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
from time import time

from geoNet.putils import get_stat_acc_plot, get_stat_vel_plot
from geoNet.utils import get_sorted_stats_code,  read_statsll

init_time = time()
#parent_dir_loc="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017/Vol1/data"
parent_dir_loc="/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105/Vol1/data"
plot_dir_accBB="accBB"
plot_dir_velBB="velBB"
#plot_dir_accLF_0pt25="accLF_0pt25"
#plot_dir_velLF_0pt25="velLF_0pt25"
#plot_dir="velBB"
loc_accBB="/".join([parent_dir_loc, plot_dir_accBB])
loc_velBB="/".join([parent_dir_loc, plot_dir_velBB])
#loc_accLF_0pt25="/".join([parent_dir_loc, plot_dir_accLF_0pt25])
#loc_velLF_0pt25="/".join([parent_dir_loc, plot_dir_velLF_0pt25])

loc_acc = loc_accBB
loc_vel = loc_velBB

#import sys
#sys.exit()
#To plot all stations in file all_stats use code below
#with open("all_stats.txt",'r') as f:
#    event_stats = f.readlines()

stats_dict = read_statsll("/nesi/projects/nesi00213/RealTime/Obs/Mw6pt6_2013-08-16_023105",
                          'geonet_stations_20170117.ll')
#from get_processed_stats_list import get_processed_stats_list
#stats_dict=get_processed_stats_list(loc_velBB,stats_dict)
#stations are sorted according to PSA
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
print("Done in {:.1f}".format(final_time - init_time))

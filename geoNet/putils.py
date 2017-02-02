"""
Convenience functions used for plotting and visualization
"""
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

from pyqc.rspectra import Response_Spectra
from pyqc.utils import get_extremum, get_stat_data, get_PSA
from pyqc.utils import get_adjusted_stat_data

def get_stat_acc_plot(loc, stat_code):

    #stat=get_stat_data(loc,stat_code)
    stat=get_adjusted_stat_data(loc,stat_code)

    fig, ax = plt.subplots(3,1)
    ext_000, argext_000 = get_extremum(stat["000"])
    ax[0].plot(stat["t"], stat["000"],c='k',label='000')
    ax[0].legend(loc='best',framealpha=1.)
    #ax[0].set_xlim(0.,250)
    ax[0].plot(stat["t"][argext_000],ext_000,marker='o')
    ax[0].annotate('{:8.4f}'.format(ext_000), xy=(stat["t"][argext_000], stat["000"][argext_000]))

    ax[1].plot(stat["t"], stat["090"],c='b',label='090')
    ax[1].legend(loc='best',framealpha=1.)
    ax[1].set_ylabel("Acceleration, a[g]")
    #ax[1].set_xlim(0.,250)
    ext_090, argext_090 = get_extremum(stat["090"])
    ax[1].plot(stat["t"][argext_090],ext_090,marker='o')
    ax[1].annotate('{:8.4f}'.format(ext_090),xy=(stat["t"][argext_090], stat["090"][argext_090]))

    ax[2].plot(stat["t"], stat["ver"],c='r',label='ver')
    ax[2].legend(loc='best',framealpha=1.)
    #ax[2].set_xlim(0.,250)
    ext, argext = get_extremum(stat["ver"])
    ax[2].plot(stat["t"][argext],ext,marker='o')
    ax[2].annotate('{:8.4f}'.format(ext), xy=(stat["t"][argext], stat["ver"][argext]))

    fig.suptitle(stat_code + ", PGA geometric mean=" + "{:.4f}".format(np.sqrt(np.abs(ext_000*ext_090))))
    return stat, fig, ax


def get_stat_vel_plot(loc, stat_code):

    #stat=get_stat_data(loc,stat_code)
    stat=get_adjusted_stat_data(loc,stat_code)
    
    fig, ax = plt.subplots(3,1)
    ax[0].plot(stat["t"], stat["000"],c='k',label='000')
    ax[0].legend(loc='best',framealpha=1.)
    #ax[0].set_xlim(0.,250)
    ext_000, argext_000 = get_extremum(stat["000"])
    ax[0].plot(stat["t"][argext_000],ext_000,marker='o')
    ax[0].annotate('{:8.2f}'.format(ext_000), xy=(stat["t"][argext_000], stat["000"][argext_000]))

    ax[1].plot(stat["t"], stat["090"],c='b',label='090')
    ax[1].legend(loc='best',framealpha=1.)
    ax[1].set_ylabel("Velocity (cm/s)")
    #ax[1].set_xlim(0.,250)
    ext_090, argext_090 = get_extremum(stat["090"])
    ax[1].plot(stat["t"][argext_090], ext_090,marker='o')
    ax[1].annotate('{:8.2f}'.format(ext_090), xy=(stat["t"][argext_090], stat["090"][argext_090]))

    ax[2].plot(stat["t"], stat["ver"],c='r',label='ver')
    ax[2].legend(loc='best',framealpha=1.)
    #ax[2].set_xlim(0.,250)
    ext, argext = get_extremum(stat["ver"])
    ax[2].plot(stat["t"][argext], ext, marker="o")
    ax[2].annotate('{:8.2f}'.format(ext), xy=(stat["t"][argext], stat["ver"][argext]))

    fig.suptitle(stat_code + ", PGV geometric mean=" + "{:.2f}".format(np.sqrt(np.abs(ext_000*ext_090))))
    return stat, fig, ax


def get_stat_PSA_plot(loc, stat_code, period=np.logspace(np.log10(0.01), np.log10(10.), 100)):
    """
    """
    #stat=get_stat_data(loc,stat_code)
    stat=get_adjusted_stat_data(loc,stat_code)
    dt = stat["t"][1] - stat["t"][0]
    PSA=get_PSA(stat, dt, period) 

    fig, ax = plt.subplots()
    ax.plot(period, PSA["000"],c='k', ls='--',label='000')
    ax.plot(period, PSA["090"],c='b', ls='-.',label='090')
    ax.plot(period, PSA["ver"],c='r', ls='dotted',label='ver')
    ax.plot(period, PSA["geom"],c='g',lw='3', ls='solid',label='geom')
    ax.legend(loc='best')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(stat_code)
    ax.set_xlabel("Period, T[s]")
    ax.set_ylabel("Spectral acc, Sa[g]")
    return stat, PSA, fig, ax

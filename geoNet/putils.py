"""
Convenience functions used for plotting and visualization
"""
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.backends.backend_pdf import PdfPages

from geoNet.rspectra import Response_Spectra
from geoNet.utils import get_extremum, get_stat_data, get_PSA
from geoNet.utils import get_adjusted_stat_data

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



def plot_bias(periods, bias, std, savefig=True, figName="bias", ext="png"):
    """
    bias, std are those returned by get_bias
    """
    fig_bias, ax_bias = plt.subplots()
    ax_bias.semilogx(periods, bias)
    ax_bias.set_xlim([0.01, 10])
    ax_bias.set_ylim(-2.5, 2.5)
    ax_bias.set_xlabel('Vibration period, T (s)')
    ax_bias.set_ylabel('ln(obs/sim)-pSA')
    plt.setp(ax_bias.get_lines(), c='k', lw='5',ls="solid")
    ax_bias.grid(b=True, axis='y',which='major')        #turn the major grid lines on for y axis
    ax_bias.grid(b=True, axis='x',which='minor')
    plt.fill_between(periods, bias-std, bias+std,
                     facecolor=[0.8,0.8,0.8], edgecolor='k',linestyle='dashed',
                     linewidth=.5, alpha=0.5)
    fig_bias.set_tight_layout(True)
    if savefig:
        fig_bias.savefig(figName+"."+ext)

    return fig_bias, ax_bias

def plot_IMvsRrup(Rrup, IM, IM_std, fig=None, ax=None):
    """
    Plots gmpe IM vs Rrup.
    Rrup:
        1d numpy array.
    IM:
        1d numpy array.
    NOTE:
        No checks are performed to verify that fig, ax are indeed matplotlib objects
    return:
        matplotlib fig, ax. If fig and ax are not None, IMvsRrup plot will be
        overlayed.
    """
    if fig is None and ax is None:
        fig, ax = plt.subplots()
    else:
        pass

    #zorder places the plot beneath all else
    ax.plot(Rrup, IM,c='k', ls='solid', zorder=-1, label='GMPE')
    ax.plot(Rrup, IM*np.exp(IM_std), c='k',ls='dashed', zorder=-1)
    ax.plot(Rrup, IM*np.exp(-IM_std), c='k', ls='dashed', zorder=-1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(Rrup[0],Rrup[-1])

    return fig, ax

def plot_SMS_IM(Rrup, IM_sim, IM_obs):
    """
    Rrup:
        1d numpy array
    IM_sim:
        1d numpy array, e.g PGV, PGA, pSA for simulations
    IM_obs:
        1d numpy array, e.g PGV, PGA, pSA for observations
    """
    fig, ax = plt.subplots()
    ax.scatter(Rrup, IM_obs,  c='g',linewidths=2.,  marker='+', s=40, label='Obs')
    ax.scatter(Rrup, IM_sim,  c='r', marker='o', s=40, label='Sim')
    ax.set_xscale('log')
    ax.set_yscale('log')
    return fig, ax

def plot_SMS_IM_ratio(Rrup, IM_ratio):
    """
    Rrup: 
        1d numpy array
    IM_ratio:
        1d numpy array 
    """
    IM_ratio = np.log(IM_ratio)
    fig, ax = plt.subplots()
    ax.semilogx(Rrup, IM_ratio, linestyle='None', linewidth=5, color='b', 
                marker='s', markersize=12)
    ax.grid(b=True, axis='y',which='major')
    ax.grid(b=True, axis='x',which='minor')
    #Note mean being passed on as median.
    biasMu=np.mean(IM_ratio)
    biasSigma=np.std(IM_ratio, ddof=1)
    biasString='median={:.2f}, sigma={:.2f}'.format(biasMu,biasSigma)
    ax.annotate(biasString,xy=(0.95,0.95),xycoords='axes fraction',
                verticalalignment='center', horizontalalignment='right', fontsize=12)
    return fig, ax


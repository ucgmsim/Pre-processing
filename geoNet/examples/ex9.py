"""
Does all the plotting done in post processing for observations only
"""
from event import keyValueFromTxt

def plot_all(event_info_fname):
    """
    """
    import numpy as np
    import matplotlib
    # Force matplotlib to not use any Xwindows backend.
    matplotlib.use('Agg')
    from matplotlib import pylab as plt
    from geoNet import utils, putils
    from geoNet.gmpe import readStationFile as rsf
    from geoNet.gmpe.calculateGMPE import set_faultprop

    info = keyValueFromTxt(event_info_fname)

    stats_dict_obs = utils.read_statsll(info["loc_statsll"], info["fname_statsll"])
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict_obs = utils.get_processed_stats_list(
                     "/".join([info["loc_V1A"], info["obs_velDir"]]),
                      stats_dict_obs, verbose=True)

    FiniteFault = rsf.readSrfFile("/".join([info["srf_dir"], info["srf_file"]]))
    FiniteFault = rsf.Points_np(FiniteFault)
    #set fault properties
    faultprop = set_faultprop(Mw=float(info["Mw"]), rake=float(info["rake"]),
                              dip=float(info["dip"]), Ztor=float(info["dip"]))

    sorted_stats_code, sms_rrup_rjb = utils.get_SMS_Rrups(stats_dict_obs, FiniteFault)
    sms_rrup, sms_rjb = sms_rrup_rjb[:,0], sms_rrup_rjb[:,1]
    #**************************************************************************
    # (3) calculate pSA for simulations and observation for the above SMSs
    #Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
    #one or the other
    g=981. #cm/s^2
    periods = np.logspace(start=np.log10(0.01), stop=np.log10(10.), num=100, base=10)

    pSA_obs = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["loc_V1A"], info["obs_accDir"]]),
                          comp='geom')


    for i, stat_code in enumerate(sorted_stats_code):
        fig, ax = plt.subplots()
        ax.loglog(periods, pSA_obs[i,:],c='k',lw='1', ls='solid',label=stat_code+' obs')
        ax.legend(loc='best')
        #ax.set_title(stat_code)
        ax.set_xlabel("Period, T[s]")
        ax.set_ylabel("Spectral acc, Sa[g]")    
        fig.savefig("figs_obs/pSA_{:s}".format(stat_code)+".png")
        fig.clear()
        plt.close('all')
    #*******************************************************************************
    #Above calculations are repeated for selected periods, this avoids having 
    #to interpolate
    ## (3) calculate pSA for simulations and observation for the above SMSs
    ##Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
    ##one or the other
    #g=981. #cm/s^2
    periods=np.array([0.1, 0.2, 0.5, 1.0, 3.0, 5.0, 8.0, 10.0])


    pSA_obs = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["loc_V1A"], info["obs_accDir"]]),
                          comp='geom')

    # (4) calculate pSA for GMPE
    Rrups_gmpe = np.logspace(np.log10(5.),np.log10(100.),30)
    pSA_gmpe, pSA_gmpe_std = utils.get_empIM_v2(Rrups_gmpe, periods, faultprop)

    # (5) (a) Plot obs, sim pSAs (b) underlay gmpe pSA
    #def plot_SMS_IM(Rrup, IM_sim, IM_obs):
    for i, T in enumerate(periods):
        fig, ax = plt.subplots()
        ax.scatter(sms_rrup, pSA_obs[:,i],  c='g',linewidths=2.,  marker='+', s=40, label='Obs')
        ax.set_xscale('log')
        ax.set_yscale('log')
        #Now underlay gmpe predictions
        fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, pSA_gmpe[:,i], pSA_gmpe_std[:,i], fig=fig, ax=ax)
        ax.legend(loc="best", scatterpoints=1)
        fig.savefig("figs_obs/pSA{:.1f}".format(T)+".png")
        fig.clear()
        plt.close('all')


    #***************************************************************************
    # (3) calculate PGV, PGA for simulations and observation for the above SMSs
    #Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
    #one or the other
    #g=981. #cm/s^2

    PGV_obs = utils.get_SMS_PGV(sorted_stats_code, 
                          "/".join([info["loc_V1A"], info["obs_velDir"]]),
                          absVal=False)

    PGA_obs = utils.get_SMS_PGA(sorted_stats_code,
                          "/".join([info["loc_V1A"], info["obs_accDir"]]),
                          absVal=False)

    # (4) calculate PGV, PGA with GMPE
    #get_empIM_v2(Rrup, period, faultprop, Rjb=None, Rtvz=0., V30measured=0., V30=250.):
    Rrups_gmpe = np.logspace(np.log10(5.),np.log10(100.),30)
    period=-1
    PGV_gmpe, PGV_gmpe_std = utils.get_empIM_v2(Rrups_gmpe, period, faultprop)
    period=0
    PGA_gmpe, PGA_gmpe_std = utils.get_empIM_v2(Rrups_gmpe, period, faultprop)
    # (5)  PGV plots
    #Note -1 is the geometric mean component
    fig, ax = plt.subplots()
    ax.scatter(sms_rrup, PGV_obs[:,-1],  c='g',linewidths=2.,  marker='+', s=40, label='Obs')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #Now underlay gmpe predictions
    fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, PGV_gmpe[:,0], PGV_gmpe_std[:,0], fig=fig, ax=ax)
    ax.legend(loc="best", scatterpoints=1)
    fig.savefig("figs_obs/PGV.png")
    plt.close('all')


    # (6)  PGA plots
    #Note -1 is the geometric mean component
    fig, ax = plt.subplots()
    ax.scatter(sms_rrup, PGA_obs[:,-1],  c='g',linewidths=2.,  marker='+', s=40, label='Obs')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #Now underlay gmpe predictions
    fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, PGA_gmpe[:,0], PGA_gmpe_std[:,0], fig=fig, ax=ax)
    ax.legend(loc="best", scatterpoints=1)
    fig.savefig("figs_obs/PGA.png")
    plt.close('all')

    return

if __name__ == "__main__":
    import argparse                                                                 
    parser = argparse.ArgumentParser(description='Use for event analyses.')
    parser.add_argument('-f','--fname', help='name of input file is required', required=True)         
    args = vars(parser.parse_args())
    plot_all(args['fname'])

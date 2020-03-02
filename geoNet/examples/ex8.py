"""
Does all the plotting done in sim & obs post processing
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
    from geoNet.geoNet import putils
    from geoNet.geoNet import utils
    from geoNet.gmpe import readStationFile as rsf
    from geoNet.gmpe.calculateGMPE import set_faultprop

    info = keyValueFromTxt(event_info_fname)

    stats_dict_obs = utils.read_statsll(info["loc_statsll"], info["fname_statsll"])
    stats_dict_sim = utils.read_statsll(info["loc_statsll_sim"], info["fname_statsll_sim"])
    #Some of the SMSs may not be processed. Assign stats_dict to only those that were processed.
    stats_dict_obs = utils.get_processed_stats_list(
                     "/".join([info["loc_V1A"], info["obs_velDir"]]),
                      stats_dict_obs, verbose=True)

    stats_dict_simObs = {key:stats_dict_sim[key] for key in (set(stats_dict_sim) &
                                                            set(stats_dict_obs))
                        }
    FiniteFault = rsf.readSrfFile("/".join([info["srf_dir"], info["srf_file"]]))
    FiniteFault = rsf.Points_np(FiniteFault)
    #set fault properties
    faultprop = set_faultprop(Mw=float(info["Mw"]), rake=float(info["rake"]),
                              dip=float(info["dip"]), Ztor=float(info["dip"]))

    sorted_stats_code, sms_rrup_rjb = utils.get_SMS_Rrups(stats_dict_simObs, FiniteFault)
    sms_rrup, sms_rjb = sms_rrup_rjb[:,0], sms_rrup_rjb[:,1]
    #**************************************************************************
    # (3) calculate pSA for simulations and observation for the above SMSs
    #Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
    #one or the other
    g=981. #cm/s^2
    periods = np.logspace(start=np.log10(0.01), stop=np.log10(10.), num=100, base=10)

    pSA_sim = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["bb_dir"], info["sim_accDir"]]),
                                comp='geom') / g
    pSA_obs = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["loc_V1A"], info["obs_accDir"]]),
                                comp='geom')
    bias, std = utils.get_bias(pSA_obs, pSA_sim, rescale=False)

    # (4) Finally plot the bias
    figBias, axBias = putils.plot_bias(periods, bias, std,
                                       savefig=True, figName="figs/bias", ext="png")

    for i, stat_code in enumerate(sorted_stats_code):
        fig, ax = plt.subplots()
        #ax.plot(period, PSA["000"],c='k', ls='--',label='000')
        #ax.plot(period, PSA["090"],c='b', ls='-.',label='090')
        #ax.plot(period, PSA["ver"],c='r', ls='dotted',label='ver')
        ax.loglog(periods, pSA_sim[i,:],c='r',lw='1', ls='solid',label=stat_code + '\nsim')
        ax.loglog(periods, pSA_obs[i,:],c='k',lw='1', ls='solid',label='obs')
        ax.legend(loc='best')
        #ax.set_title(stat_code)
        ax.set_xlabel("Period, T[s]")
        ax.set_ylabel("Spectral acc, Sa[g]")    
        fig.savefig("figs/pSA_{:s}".format(stat_code)+".png")
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


    pSA_sim = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["bb_dir"], info["sim_accDir"]]),
                                comp='geom') / g
    pSA_obs = utils.get_SMS_pSA(sorted_stats_code, periods,
                          "/".join([info["loc_V1A"], info["obs_accDir"]]),
                                comp='geom')

    # (4) calculate pSA for GMPE
    Rrups_gmpe = np.logspace(np.log10(5.),np.log10(100.),30)
    pSA_gmpe, pSA_gmpe_std = utils.get_empIM_v2(Rrups_gmpe, periods, faultprop)

    # (5) (a) Plot obs, sim pSAs (b) underlay gmpe pSA
    #def plot_SMS_IM(Rrup, IM_sim, IM_obs):
    for i, T in enumerate(periods):
        fig, ax = putils.plot_SMS_IM(sms_rrup, pSA_sim[:, i], pSA_obs[:, i])
        #Now underlay gmpe predictions
        fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, pSA_gmpe[:, i], pSA_gmpe_std[:, i], fig=fig, ax=ax)
        ax.legend(loc="best", scatterpoints=1)
        fig.savefig("figs/pSA{:.1f}".format(T)+".png")
        fig.clear()
        plt.close('all')

    # (4)  Plot ratio of obs to sim pSAs
    for i, T in enumerate(periods):
        fig, ax = putils.plot_SMS_IM_ratio(sms_rrup, pSA_obs[:, i] / pSA_sim[:, i])
        ax.set_ylim(-2.,2.)
        fig.savefig("figs/pSA{:.1f}_ratio".format(T)+".png")
        fig.clear()
        plt.close('all')

    #***************************************************************************
    # (3) calculate PGV, PGA for simulations and observation for the above SMSs
    #Note: Observations acc has g units but simulations acc is in cm/s^2. Rescale
    #one or the other
    #g=981. #cm/s^2

    PGV_sim = utils.get_SMS_PGV(sorted_stats_code,
                          "/".join([info["bb_dir"], info["sim_velDir"]]),
                                absVal=False)
    PGV_obs = utils.get_SMS_PGV(sorted_stats_code, 
                          "/".join([info["loc_V1A"], info["obs_velDir"]]),
                                absVal=False)

    PGA_sim = utils.get_SMS_PGA(sorted_stats_code,
                          "/".join([info["bb_dir"], info["sim_accDir"]]),
                                absVal=False) / g
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
    fig, ax = putils.plot_SMS_IM(sms_rrup, PGV_sim[:, -1], PGV_obs[:, -1])
    #Now underlay gmpe predictions
    fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, PGV_gmpe[:, 0], PGV_gmpe_std[:, 0], fig=fig, ax=ax)
    ax.legend(loc="best", scatterpoints=1)
    fig.savefig("figs/PGV.png")
    plt.close('all')

    fig, ax = putils.plot_SMS_IM_ratio(sms_rrup, PGV_obs[:, -1] / PGV_sim[:, -1])
    ax.set_ylim(-2.,2.)
    fig.savefig("figs/PGV_ratio.png")
    plt.close('all')


    # (6)  PGA plots
    #Note -1 is the geometric mean component
    fig, ax = putils.plot_SMS_IM(sms_rrup, PGA_sim[:, -1], PGA_obs[:, -1])
    #Now underlay gmpe predictions
    fig, ax = putils.plot_IMvsRrup(Rrups_gmpe, PGA_gmpe[:, 0], PGA_gmpe_std[:, 0], fig=fig, ax=ax)
    ax.legend(loc="best", scatterpoints=1)
    fig.savefig("figs/PGA.png")
    plt.close('all')

    fig, ax = putils.plot_SMS_IM_ratio(sms_rrup, PGA_obs[:, -1] / PGA_sim[:, -1])
    ax.set_ylim(-2.,2.)
    fig.savefig("figs/PGA_ratio.png")
    plt.close('all')
    return

if __name__ == "__main__":
    import argparse                                                                 
    parser = argparse.ArgumentParser(description='Use for event analyses.')
    parser.add_argument('-f','--fname', help='name of input file is required', required=True)         
    args = vars(parser.parse_args())
    plot_all(args['fname'])

"""
Convenience functions used for reading, writting strong motion station 
data, station list files (.ll) and more.
"""
import numpy as np
import os
from glob import glob
import datetime
#https://docs.python.org/2.5/whatsnew/pep-328.html
#from quakecore.rspectra import Response_Spectra
#from quakecore.gmpe.Bradley_2010_Sa import Bradley_2010_Sa


def readGP(loc, fname):
    """
    Convinience function for reading files in the Graves and Pitarka format
    """
    with open("/".join([loc, fname]), 'r') as f:
        lines = f.readlines()
    
    data = []

    for line in lines[2:]:
        data.append([float(val) for val in line.split()])

    data=np.concatenate(data) 
    
    line1=lines[1].split()
    num_pts=float(line1[0])
    dt=float(line1[1])
    shift=float(line1[4])

    return data, num_pts, dt, shift


def get_GP_header(stat_code, size, delta_t, time_delay):
    """
    Return header for GP file
    """

    header_000 = stat_code + " 0 broadband\n"
    header_090 = stat_code + " 90 broadband\n"
    header_ver = stat_code + " ver broadband\n"
    stat_info= ("{:<10d}"+ 7*"{:<10.3f}"+"\n").format(size, delta_t, 0., 0. ,time_delay, 0., 0., 0.)

    header_000+=stat_info
    header_090+=stat_info
    header_ver+=stat_info

    return header_000, header_090, header_ver


def writeGP(loc, fname, data, header, ncol=6):
    """
    Convinience function for writing files in the Graves and Pitarka format
    """
    size = len(data)
    nrow = size / ncol
    size_last_row = size % ncol
    
    lines = ""
    for line in np.reshape(xrange(nrow*ncol), (nrow, ncol)):
        for val in line:
            lines += "{:^20.6e}".format(data[val]) + 3*" "
        lines = lines.rstrip(3*" ") + "\n"

    if size_last_row:
        for i in xrange(1, size_last_row+1):
            lines += "{:^20.6e}".format(data[-i]) + 3*" "
    lines = lines.rstrip(3*" ")

    with open("/".join([loc, fname]), 'w') as f:
        f.writelines(header)
        f.writelines(lines)
    return


def writeGP_stat_data(loc, stat_data, time_delay, ncol=6):
    """
    writes station data:  000, 090, ver components
    stat_data:
        is of type returned by get_stat_data
    """
    size=stat_data['000'].size
    delta_t=stat_data['t'][1]-stat_data['t'][0]
    stat_code=stat_data['name']
    header_000, header_090, header_ver = get_GP_header(
                                       stat_code, size, delta_t, time_delay)       

    writeGP(loc, stat_code+".000", stat_data['000'], header_000, ncol)
    writeGP(loc, stat_code+".090", stat_data['090'], header_090, ncol)
    writeGP(loc, stat_code+".ver", stat_data['ver'], header_ver, ncol)

    return

def adjust_for_time_delay(ts, dt, shift):
    """
        ts: time series data
    """
    t0_index = int(shift/dt)
    if t0_index == 0:
        num_pts = ts.size
    elif t0_index > 0:
        ts = np.concatenate((np.zeros(t0_index), ts))
        num_pts = ts.size
    elif t0_index <0:
        ts = ts[np.abs(t0_index):]
        num_pts = ts.size

    return ts, num_pts, dt


def get_adjusted_stat_data(loc, stat_code):
    """
    Time series data is adjust for time shift

    returns a dictionary with components:
        000, 090, ver, t, stat_code
    """
    stat_data = {"000": None, "090": None, "ver":None, "t": None, "name": None}
    g=981. #cm/s^2
    stat_data["000"], num_pts, dt, shift  = readGP(loc, ".".join([stat_code, "000"]))
    stat_data["090"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "090"]))   
    stat_data["ver"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "ver"]))


    stat_data["000"], num_pts, dt = adjust_for_time_delay(stat_data["000"], dt, shift)
    stat_data["090"], num_pts, dt = adjust_for_time_delay(stat_data["090"], dt, shift)
    stat_data["ver"], num_pts, dt = adjust_for_time_delay(stat_data["ver"], dt, shift)

    t = np.arange(num_pts)*dt
    stat_data["t"] = t

    stat_data["name"]=stat_code
    return stat_data


def get_stat_data(loc, stat_code):
    """
    returns a dictionary with components:
        000, 090, ver, t, stat_code
    """
    stat_data = {"000": None, "090": None, "ver":None, "t": None, "name": None}
    g=981. #cm/s^2
    stat_data["000"], num_pts, dt, shift  = readGP(loc, ".".join([stat_code, "000"]))
    stat_data["090"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "090"]))   
    stat_data["ver"], num_pts, dt, shift = readGP(loc, ".".join([stat_code, "ver"]))


    #stat_data["090"] /=g
    #stat_data["ver"] /=g

    t = np.arange(num_pts)*dt
    stat_data["t"] = t

    stat_data["name"]=stat_code
    return stat_data


def get_event_data(loc, stats_dict):
    """
    loc: location of where data files are
    """
    for stat_code in stats_dict.keys():
        yield get_stat_data(loc, stat_code)


def get_extremum(x):
    """
    x must be numpy array
    """
    min_val = x.min()
    argmin = x.argmin()
    max_val = x.max()
    argmax = x.argmax()
    extremum = min_val
    argextremum=argmin
    if abs(min_val) < abs(max_val):
        extremum = max_val
        argextremum=argmax

    return extremum, argextremum


def get_PSA(acc_dict, dt, periods,
                 xi=0.05,  m=1, gamma=0.5, beta=0.25):
    """
    acc is a directionary of acceleration time series of all three components
    returns:
        dictionary pSA
    """

    try:
        pSA_000 = Response_Spectra(acc_dict["000"], dt, xi, periods, m, gamma, beta)
        pSA_090 = Response_Spectra(acc_dict["090"], dt, xi, periods, m, gamma, beta)
        pSA_ver = Response_Spectra(acc_dict["ver"], dt, xi, periods, m, gamma, beta)
        pSA = {"000":pSA_000, "090":pSA_090, "ver":pSA_ver, "geom":np.sqrt(pSA_000*pSA_090)}
    except Exception as e:
        print(e)

    return  pSA


def get_event_PSA(event_acc, period, xi=0.05, m=1., gamma=0.5, beta=0.25):
    """
    event_acc: an iterator of event acceleration data
    """
    for acc_data in event_acc:
        dt = acc_data["t"][1] - acc_data["t"][0]
        PSA=get_PSA(acc_data, dt, period, xi=xi, m=m, gamma=gamma, beta=beta) 
        PSA["name"]=acc_data["name"]
        yield PSA



def read_statsll(loc, name):
    """ 
    fname:
        .ll format station file that contains lon, lat, stat_code
    returns a dictionary
    """
    #lon, lat, stat_code = np.genfromtxt(fname, dtype="%10.4f, %10.4f, %5s", unpack=True)
    #stats_dict = {"lon":lon, "lat":lat, "stat":stat_code}

    stats_dict = {}
    with open("/".join([loc, name]), 'r') as f:
        for line in f:
            lon, lat, stat_code = line.split()
            lon = float(lon)
            lat = float(lat)
            stat_code = stat_code.strip() #remove space and newline character
            stats_dict[stat_code] = (lon, lat)

    return stats_dict


def get_sorted_stats_code(loc, stats_dict, comp="geom"):
    """
        Applies for vel and acc data only
        loc: location of where .000, .090, .ver data files are
        returns: 
            A list. Each list item is a dictionary
    """
    sorted_stats_code = []
    
    for stat_code in stats_dict.keys():
        stat_code = stat_code.strip() #remove spaces and \n
        stat_data = get_stat_data(loc, stat_code)
        max_000 = np.max(np.abs(stat_data["000"]))
        max_090 = np.max(np.abs(stat_data["090"]))
        max_geom= np.sqrt(max_000*max_090)
        sorted_stats_code.append({"000":max_000, "090":max_090, "geom":max_geom, "name":stat_code})

    from operator import itemgetter
    sorted_stats_code = sorted(sorted_stats_code, key=itemgetter('geom'))
    return sorted_stats_code 


def get_emp_pSA(periods, siteprop, faultprop):
    """
    Note that siteprop period is changed
    returns:
        PSA_emp where np.shape(PSA_emp) == (len(periods), 2)
    """

    if np.isscalar(periods):
        periods = [periods]

    periods = np.asarray(periods)
    pSA_emp = np.empty((periods.size, 2))

    for i, T in enumerate(periods):
        siteprop.period=T 
        pSA_emp[i,0], pSA_emp_sigm = Bradley_2010_Sa(siteprop,faultprop)
        pSA_emp[i,1] = pSA_emp_sigm[0]

    return pSA_emp


def get_events_stations(fname_all_geoNet_stats="all_geoNet_stats.ll",
                        loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo",
                        loc_Vol1="/".join([os.getcwd(), "Vol1"]),
                        save_stats=False, fname=None, loc=os.getcwd()):

    all_geoNet_stats = read_statsll(loc_all_geoNet_stats, fname_all_geoNet_stats)

    event_stats_V1A = glob("/".join([loc_Vol1,"data","*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    
    event_stats= {}
    for V1A_file in event_stats_V1A:
        year, event_id, stat_code, V1A = V1A_file.split("_")
        event_stats[stat_code] = (None, None)
        if all_geoNet_stats.has_key(stat_code):
            event_stats[stat_code] = all_geoNet_stats[stat_code]

    if save_stats == True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        if fname is None:
            fname = "_".join([year, event_id, "eventStats", str(datetime.date.today())])
            fname+=".ll"
            with open("/".join([loc,fname]),'w') as f:
                for key, value in event_stats.items():
                    if value[0] is None:
                        print("{:10s} not found in all_geoNet_stats, add this manually to event_stats.ll".format(key))
                    else:
                        line="{:10.4f} {:10.4f} {:10s}".format(value[0], value[1], key)
                    f.write(line+"\n")

    return event_stats

def vel2acc(timeseries, dt):
    """
    Source is shared_ts.py
    Differentiate following Rob Graves' code logic.
    """
    return np.diff(np.hstack(([0], timeseries)) * (1.0 / dt))

def acc2vel(timeseries, dt):
    """
    Source is shared_ts.py
    Integrates following Rob Graves' code logic (simple).
    """
    return np.cumsum(timeseries) * dt

"""
Convenience functions used for reading, writting strong motion station 
data, station list files (.ll) and more.
"""
import numpy as np
import os
from glob import glob
import datetime
from itertools import izip
from scipy.interpolate import UnivariateSpline as US
from scipy.integrate import cumtrapz
from scipy import signal, fftpack
#https://docs.python.org/2.5/whatsnew/pep-328.html
from geoNet.rspectra import Response_Spectra
from geoNet.gmpe.Bradley_2010_Sa import Bradley_2010_Sa
from geoNet.gmpe.calculateGMPE import set_siteprop
from geoNet.gmpe import readStationFile as rsf

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


def get_GP_header(stat_code, size, delta_t, time_delay, comment="broadband"):
    """
    Return header for GP file
    """

    #header_000 = stat_code + " 0 broadband\n"
    #header_090 = stat_code + " 90 broadband\n"
    #header_ver = stat_code + " ver broadband\n"
    comment = str(comment)
    comment = comment.strip("\n")
    header_000 = "{:s} 000 {:s}\n".format(stat_code, comment)
    header_090 = "{:s} 090 {:s}\n".format(stat_code, comment)
    header_ver = "{:s} ver {:s}\n".format(stat_code, comment)
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
    #for stat_code in stats_dict.keys():
    for stat_code in stats_dict:
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

def write_statsll(loc, fname, stats_dict):
    """
    loc:
        where fname is saved
    stats_dict:
        is of type returned by read_statsll
    """

    with open("/".join([loc,fname]),'w') as f:
        for stat_code in stats_dict:
            (lon, lat) = stats_dict[stat_code]
            f.write("{:<15.4f} {:^15.4f} {:^10s} \n".format(lon, lat, stat_code))

    return
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


def get_empIM(periods, siteprop, faultprop):
    """ 
    This is a wrapper function to Bradley_2010_Sa. For 
    a fixed Rrup it calls Bradley_2010_Sa periods.size times.
    Note:
        In siteprop period is changed. Thus when setting siteprop one can 
        set siteprop.period=np.nan
    returns:
        PGV, PGA, pSA with np.shape(IM) = (periods.size,2)
    """

    if np.isscalar(periods):
        periods = [periods]

    periods = np.asarray(periods)
    empIM = np.empty((periods.size, 2)) 

    for i, T in enumerate(periods):
        siteprop.period=T 
        empIM[i,0], empIM_sigm = Bradley_2010_Sa(siteprop,faultprop)
        empIM[i,1] = empIM_sigm[0]

    return empIM


def get_empIM_v2(Rrup, period, faultprop, Rjb=None, Rtvz=0., V30measured=0., V30=250.):
    """
    Similar to get_empIM except that both Rrup, period and V30 can be a numpy arrays or scalar.
    Note that in set_siteprop we have 
        self.Rjb=np.sqrt(np.max((0,Rrup**2-faultprop.Ztor**2)))
    and therefore Rjb has no effect currently
    """
    if np.isscalar(Rrup):
        Rrup = [Rrup]
    if np.isscalar(period):
        period = [period]

    if np.isscalar(V30):
        V30 = [V30]*len(Rrup)

    Rrup=np.asarray(Rrup)
    if Rjb is None: 
        Rjb=-Rrup
    V30=np.asarray(V30)
    if not Rrup.size == V30.size: raise AssertionError("Rrup.size not equal to V30.size")
    period=np.asarray(period)


    empIM_values = np.empty((Rrup.size, period.size))
    empIM_sigmas = np.empty((Rrup.size, period.size))
    for i, (x,y,z) in enumerate(izip(Rrup, Rjb, V30)):
        #x,y = Rrup[i], Rjb[i]
        siteprop = set_siteprop(faultprop, period=np.nan, Rrup=x, Rjb=y,
                                Rtvz=Rtvz, V30measured=V30measured, V30=z)

        empIM = get_empIM(period, siteprop, faultprop)
        empIM_values[i] = empIM[:,0]
        empIM_sigmas[i] = empIM[:,1]

    return empIM_values, empIM_sigmas


def get_processed_stats_list(loc, stats_dict, verbose=False):
    """
    loc:
        location of .000, .090, .ver files
    stats_dict:
        Dictionary of stats as read by read_statsll function
    return:
        Dictionary of stats found at loc
    """

    processed_stats_codes = glob("/".join([loc,"*.000"]))
    psc=processed_stats_codes

    #remove path name, the for stat_code.ext remove ext
    psc = [os.path.basename(_).split(".")[0] for _ in psc]

    processed_stats_dict = {}
    stats_not_in_stats_dict= []
    for stat_code in psc:
        if stats_dict.has_key(stat_code):
            processed_stats_dict.update({stat_code:
                                        stats_dict[stat_code]})
        else:
            stats_not_in_stats_dict.append(stat_code)

    if len(stats_not_in_stats_dict) > 0 and verbose:
        print("(lon, lat) for the following stations were not in stats_dict:\n")
        print(stats_not_in_stats_dict)

    return processed_stats_dict


def statcode_lonlat_from_V1A(loc,fname):
    """
    loc: location of where geoNet data file (.V1A) is located
    fname: name of .V1A file
    Note: See http://info.geonet.org.nz/display/appdata/Accelerogram+Data+Filenames+and+Formats
        for geoNet file formats
    """
    with open ("/".join([loc,fname]), 'r') as f:
        for line_num, line in enumerate(f):
            if line_num == 1: line2 = line
            if line_num == 18:
                line19 = line
                break

    #assume line 2 is of the form:
    #Site CTZ       43 44 14S  176 37 03W
    stat_code = line2.split()[1]
    line19 = line19.split()
    lat = -(float(line19[0]) + float(line19[1])/60 + float(line19[2])/3600)
    lon = float(line19[3]) + float(line19[4])/60 + float(line19[5])/3600

    return {stat_code: (lon, lat)} #{stat_code: (round(lon,4), round(lat,4))}



def parse_V1A_fname(fname):
    """
    Note: See http://info.geonet.org.nz/display/appdata/Accelerogram+Data+Filenames+and+Formats
    The file name can take three forms:
    XNNSITEJJ, where X is the instrument code1, NN is the two-digit year 
    of acquisition, SITE is the 4 character site code, and JJ is the site's
    unique accelerogram number for that year.
    YYYYMMDD_HHMM_SITE, where YYYY/MM/DD HH:MM is the earthquake origin time 
    (UTC) and SITE is the 3 or 4 character site code.
    YYYYMMDD_HHMM_SITE_XX, where YYYY/MM/DD HH:MM is the instrument trigger time
    (UTC), SITE is the 3 or 4 character site code, and XX is the sensor location
    code. The latter is most appropriate for arrays of borehole sensors.

    fname: name of geoNet data file, must not include path
    """
    form1 = "XNNSITEJJ"
    form2 = "YYYYMMDD_HHMM_SITE"
    form3 = "YYYYMMDD_HHMM_SITE_XX"

    #Remove .V1A, .V2A extension
    fname = fname.split(".")[0]
    YYYYMMDD = ""
    HHMM = ""
    if fname.count("_") == form1.count("_"):
        stat_code = fname[3:7]
    elif fname.count("_") == form2.count("_"):
        YYYYMMDD, HHMM, stat_code = fname.split("_")
    elif fname.count("_") == form3.count("_"):
        YYYYMMDD, HHMM, stat_code, XX = fname.split("_")
    else:
        raise Exception("{:s} is unknow file name format for .V1A files\n".format(fname))

    return YYYYMMDD, HHMM, stat_code

def statsll_from_V1A(loc_V1A, save_stats=False,
                     fname=None, loc=os.getcwd()):
    """
    loc_V1A: location where .V1A geoNet data fils are located
    save_stats: If save_stats==True, save  fname.ll at loc. Note that 
    fname and loc are not used if save_stats==False
    Note:
        Here (lon, lat) are extracted from .V1A file
    return:
        (1) A dictionary containing station code, lon and lat in the format of 
        read_statsll. 
        (2) station file name

    """
    event_stats_V1A = glob("/".join([loc_V1A,"*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    event_stats = {}
    for V1A_file in event_stats_V1A:
        try:
            event_stats.update(
            statcode_lonlat_from_V1A(loc_V1A, V1A_file)
            )
        except Exception as e:
            print(e)
            print("Stations list will not include {:s}\n".format(V1A_file))
            continue

    #year = ""
    #event_id = ""
    for V1A_file in event_stats_V1A:
        try:
            #year, event_id, stat_code, V1A = V1A_file.split(".")[0].split("_")
            #split_file_name = V1A_file.split(".")[0].split("_")
            #year, event_id, stat_code = split_file_name[0:3]
            YYYYMMDD, HHMM, stat_code = parse_V1A_fname(V1A_file)
        except Exception as e:
            print(e)
            continue
        else:
            if YYYYMMDD is "":
                continue
            else:
                break

    if fname is None:
        fname = "_".join([YYYYMMDD, "eventStats", str(datetime.date.today())])
        fname+=".ll"

    if save_stats is True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        with open("/".join([loc,fname]),'w') as f:
            for stat_code, (lon, lat) in event_stats.items():
                f.write("{:^10.4f} {:^10.4f} {:^10s} \n".format(lon, lat, stat_code))

    return event_stats, fname


def statsll_from_V1A_v2(loc_all_geoNet_stats, fname_all_geoNet_stats, 
                        loc_V1A, save_stats=False,
                        fname=None, loc=os.getcwd()):
    """
    Note: This functions assigns (lon, lat) values after reading them from 
    loc_all_geoNet_stats. If save_stats is True then fname should be suppllied.
    loc_all_geoNet_stats must be in WGS84 coordinates as in 
    https://magma.geonet.org.nz/delta/app
    return:
        (1) A dictionary containing station code, lon and lat in the format of 
        read_statsll. 
    """

    all_geoNet_stats = read_statsll(loc_all_geoNet_stats, fname_all_geoNet_stats)

    event_stats_V1A = glob("/".join([loc_V1A,"*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    event_stats = {}
    stats_not_in_all_geoNet_stats=[]
    for V1A_file in event_stats_V1A:
        YYYYMMDD, HHMM, stat_code = parse_V1A_fname(V1A_file)
        if all_geoNet_stats.has_key(stat_code):
            event_stats[stat_code] = all_geoNet_stats[stat_code]
        else:
            stats_not_in_all_geoNet_stats.append(stat_code)

    if len(stats_not_in_all_geoNet_stats) > 0:
        print("(lon, lat) for the following stations were not in {:s}:\n".format(fname_all_geoNet_stats))
        print(stats_not_in_all_geoNet_stats)

    if fname is None:
        fname = "_".join([YYYYMMDD, "eventStats", str(datetime.date.today())])
        fname+=".ll"

    if save_stats is True:
        #assert fname is not None, "Specify name of station file to save"
        #assert loc is not None, "Specify location for station file to save"
        with open("/".join([loc,fname]),'w') as f:
            for stat_code, (lon, lat) in event_stats.items():
                f.write("{:^10.4f} {:^10.4f} {:^10s} \n".format(lon, lat, stat_code))

    return event_stats, fname

def get_events_stations(fname_all_geoNet_stats=None,
                        loc_all_geoNet_stats=None,
                        loc_Vol1="/".join([os.getcwd(), "Vol1"]),
                        save_stats=False, fname=None, loc=os.getcwd()):
    """
    This function should no longer be used
    """

    all_geoNet_stats = read_statsll(loc_all_geoNet_stats, fname_all_geoNet_stats)

    event_stats_V1A = glob("/".join([loc_Vol1,"data","*.V1A"]))
    event_stats_V1A = [os.path.basename(_) for _ in event_stats_V1A]
    
    event_stats= {}
    for V1A_file in event_stats_V1A:
        #year, event_id, stat_code, V1A = V1A_file.split(".")[0].split("_")
        split_file_name = V1A_file.split(".")[0].split("_")
        year, event_id, stat_code = split_file_name[0:3]
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

    return event_stats, fname




def get_SMS_Rrups(stats_dict, finiteFault):
    """
    stats_dict:
        of type read_statsll
    finiteFault:
        finiteFault = readStationFile.readSrfFile(srf_fname)
        finiteFault = readStationFile.Points_np(FiniteFault)
    return:
        sorted_stat_codes, array(Rrup, Rjb) 
    """
    stat_codes = stats_dict.keys()
    sms_rrups = np.empty((len(stat_codes), 2), dtype='float64')
    for i, stat_code in enumerate(stat_codes):
        lon, lat = stats_dict[stat_code]

        site = rsf.Points()
        site.Lon.append(lon)
        site.Lat.append(lat)
        site.Depth.append(0.)
        site = rsf.Points_np(site)

        Rrup, Rjb = rsf.computeSourcetoSiteDistance_np(finiteFault, site)
        sms_rrups[i,0]=Rrup
        sms_rrups[i,1]=Rjb

    ind = np.argsort(sms_rrups[:,0])
    sms_rrups[:,0] = sms_rrups[ind,0]
    sms_rrups[:,1] = sms_rrups[ind,1]
    #list don't recognize numpy slices
    sorted_stat_codes = [stat_codes[i] for i in ind]

    return sorted_stat_codes, sms_rrups

def get_SMS_pSA(stat_codes, period, loc, comp='geom',
                xi=0.05, m=1., gamma=0.5, beta=0.25):
    """
    Note: If stat_codes list is not sorted use get_SMS_Rrups to get a 
    list of stat_codes sorted according to Rrup values
    stat_codes:
        list containing SMS codes
    period:
        float scalar or array 
    loc:
        location of .000, .090 and .ver data files
    comp: 
        pSA component. Must be one of  ["000", "090", "ver", "geom"]
    return:
        array(len(stat_codes), len(period))

    """
    if not isinstance(stat_codes, list):
        raise TypeError("Expects stat_cdes to be a list")

    event_data = get_event_data(loc, stat_codes)

    if np.isscalar(period):
        period = [period]

    period = np.asarray(period)
    SMS_pSA = np.empty((len(stat_codes), period.size), dtype="float64")
    for i, stat_data in enumerate(event_data):
        dt = stat_data["t"][1] - stat_data["t"][0]
        pSA =get_PSA(stat_data, dt, period, 
                     xi=xi, m=m, gamma=gamma, beta=beta)
        SMS_pSA[i,:] = pSA[comp]

    return SMS_pSA

def get_SMS_PGV(stat_codes, loc, absVal=True):
    """
    Note: If stat_codes list is not sorted use get_SMS_Rrups to get a 
    list of stat_codes sorted according to Rrup values
    stat_codes:
        list containing SMS codes
    loc:
        location of .000, .090 and .ver data files
    return:
        array(len(stat_codes), 4) with columns 
        for PGV components 000, 090, ver, geom
    """
    if not isinstance(stat_codes, list):
        raise TypeError("Expects stat_cdes to be a list")

    event_data = get_event_data(loc, stat_codes)
    SMS_PGV = np.empty((len(stat_codes), 4), dtype="float64")
    for i, stat_data in enumerate(event_data):
        #ext for extremum, argext for argument/index of extremum
        ext_000, argext_000 = get_extremum(stat_data["000"])
        ext_090, argext_090 = get_extremum(stat_data["090"])
        ext_ver, argext_ver = get_extremum(stat_data["ver"])
        SMS_PGV[i,0] = ext_000
        SMS_PGV[i,1] = ext_090
        SMS_PGV[i,2] = ext_ver
        SMS_PGV[i,3] = np.sqrt(np.abs(ext_000*ext_090))
    
    if absVal:
        SMS_PGV[0:3] = np.abs(SMS_PGV[0:3])

    return SMS_PGV

def get_SMS_PGA(stat_codes, loc, absVal=True):
    return get_SMS_PGV(stat_codes, loc, absVal=absVal)

def get_bias(pSA_obs, pSA_sim, rescale=True):
    """
    pSA_obs, pSA_sim are numpy arrays returned by get_SMS_pSA
        array(len(stat_codes), len(period))
    rescale:
        If True rescales pSA_sim from cm/^2 to units in g
    """
    if rescale:
        g=981. #cm/s^2
        pSA_sim /= g
    ratio = np.log(pSA_obs/pSA_sim)
    #ratio has rows (for stats) and columns (for periods). Sum accross rows
    #bias = np.median(ratio, axis=0)
    bias = np.mean(ratio, axis=0)
    std  = np.std(ratio, axis=0, ddof=1)
    return bias, std





def int_stat_data_sp(stat_data):
    """
    stat_data is of type returned by get_stat_data
    return:
        integral obtained from a spline representation
    """
    int_stat_data = dict.fromkeys(stat_data.keys())
    int_stat_data['name'] = stat_data['name']
    int_stat_data['t'] = stat_data['t']

    for key in ["000", "090", "ver"]:
        #s=0 means no smoothing, and the spline passes through all data points
        #k=3 cubic interpolation
        spl = US(stat_data["t"], stat_data[key], s=0, k=3)
        ispl=spl.antiderivative(n=1)
        #In general antiderivative sets the initial value of the integral to zero 
        int_stat_data[key] = ispl(stat_data['t'])

    return int_stat_data





def int_stat_data(stat_data):
    """
    stat_data is of type returned by get_stat_data
    return:
        integral obtained from cumtrapz
    """
    int_stat_data = dict.fromkeys(stat_data.keys())
    int_stat_data['name'] = stat_data['name']
    int_stat_data['t'] = stat_data['t']

    for key in ["000", "090", "ver"]:
        int_stat_data[key] = cumtrapz(y=stat_data[key],
                                      dx=stat_data["t"][1]-stat_data["t"][0],
                                      initial=0.)
    return int_stat_data

def diff_stat_data(stat_data):
    """
    stat_data is of type returned by get_stat_data
    return:
        derivative is  obtained numpy.gradient  that does second order central
        difference.
    """
    diff_stat_data = dict.fromkeys(stat_data.keys())
    diff_stat_data['name'] = stat_data['name']
    diff_stat_data['t'] = stat_data['t']

    for key in ["000", "090", "ver"]:
        dt=stat_data["t"][1]-stat_data["t"][0]

        #axis keyword added to numpy version 1.11
        #diff_stat_data[key] = np.gradient(stat_data[key],
        #                                  dt, edge_order=2, axis=None)
        diff_stat_data[key] = np.gradient(stat_data[key],
                                          dt, edge_order=2)
    return diff_stat_data

def filt_stat_data(stat_data,freq, btype, output='sos', order=4,
                   worN=512):
    """
    Note:
        requires scipy version 15. or greater. Digital filters only
        plot(w, abs(h)) for frequency response.
        sosfilt may be replace with sosfiltfilt for zero phase.
    stat_data:
        is of type returned by get_stat_data. stat_data is modified inplace
    freq:
       scalar or length 2 sequence e.g [f_lowcut, f_highcut] 
    btype:
        {lowpass, highpass, bandpass, bandstop}
    return:
        sos
        b, a (numerator, denominator of filter)
        w, h (frequency and frequency response)
    """
    dt = stat_data['t'][1]-stat_data['t'][0]
    #sampling frequency fs
    fs = 1/dt
    #Nyquist frequncy Nyq
    Nyq = fs/2.
    if np.isscalar(freq):
        freq = [freq]
    freq = np.asarray(freq)
    Wn = freq/Nyq

    #sos = signal.butter(order, Wn, btype, analog=False)
    #signal.butter is just a wrapper function to iirfilter
    sos = signal.iirfilter(order, Wn, rp=None, rs=None, btype=btype,
                           analog=False, ftype='butter', output=output)
    for comp in ['000', '090', 'ver']:
        stat_data[comp] = signal.sosfiltfilt(sos, stat_data[comp])
        #stat_data[comp] = signal.sosfilt(sos, stat_data[comp])

    #get single transfer function from series of 2nd-order sections
    b, a = signal.sos2tf(sos)
    #get frequency response of the filter
    #default worN=512 frequencies between 0 and pi
    w, h = signal.freqz(b, a, worN=worN, whole=False, plot=None)
    #un-normalize by Nyq and convert to hz from radians/second (omega=2pif,f = omega/(2pi))
    #Note that 0<= w <=pi not 2pi
    w = w*Nyq/(np.pi)


    return {"sos":sos, 'b':b, 'a':a, "w":w, 'h':h}


def fft_stat_data(stat_data):
    """
    Note: the size of the fft is set by convention see
    https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.fftpack.fft.html#scipy.fftpack.fft
    returns:
        stat_data type with keys ['000', '090', 'ver', 'f', 'name']
    """
    fft_stat_data = {}
    for comp in ['000', '090', 'ver']:
        fft_stat_data[comp] = fftpack.fft(stat_data[comp])

    dt = stat_data['t'][1]-stat_data['t'][0]
    fft_stat_data['f'] = fftpack.fftfreq(stat_data['000'].size, d=dt)
    fft_stat_data['name']=stat_data['name']

    return fft_stat_data


def keyValueFromTxt(fname):
    """
    Parses file that has the form key=value and returns a dictionary
    """
    keyValue = dict()
    fname = os.path.abspath(fname)
    print("Reading input from {:s}\n".format(fname))
    with open(fname, 'r') as f:
        for line in f:
            #remove white spaces
            if line.startswith("#") or line.startswith("%"):
                continue
            if line in ["\n", "\r", "\rn"]:
                continue
            for delim in ['#', '%']:
                line = line.partition(delim)[0]
            line = line.strip()
            line = line.replace(" ", "")
            line = line.replace("\"", "")
            key, value = line.split("=")
            keyValue[key] = value

    return keyValue


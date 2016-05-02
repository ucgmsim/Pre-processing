_sample_file_header_format=\
"""\
Uncorrected accelerogram 20160214_001343_NBLC_20 GNS Science
Site NBLC      43 30 25S  172 43 53E     Cusp file:  20160214_001345_NBLC
New Brighton Library
Instrument Cusp_323
Resolution: 18-bit   Instr Period: 0.0050 s   Damping: 0.700
Accelerogram 20160214_001343_NBLC_20                 Processed 2016 February 22
11 km east of CHRISTCHURCH
2016 February 14         00:13:43 UT
Epicentre  43 29 50S 172 45 18E  Bearing N61E  Dist    2km  Depth    8km M  5.70
Number of points  26385    Duration 131.92 s
Raw readings in units of 1.00 mm/s/s  at intervals of 0.005 s
Data is unfiltered
Component S14E
Acceleration:  Peak  -2080.5 mm/s/s at  23.28 s     RMS   145.65 mm/s/s
Velocity record unevaluated
Displacement record unevaluated
    2016       2      14       0      13     430       0       0    2016       2
      43      29      50     172      45      18       8       0      14       0
      43      30      25     172      43      53     166     166      61       2
   26385       0       0   26385       0       0       0       0      13   25845
   200.0    0.70    0.00   0.000  0.0000   0.000   0.000   1.000     0.0    0.00
  43.507 172.731  43.497 172.755    5.70    0.00    0.00    0.00      0.      0.
  131.92    0.00    0.00  131.92  0.0050  0.0050      0.      0.      0.  9806.6
 -2080.5   23.28   145.7      0.      0.     0.0    0.00     0.0     0.0   0.000
    0.00    0.00    0.00    0.00   0.000    0.00    0.00    0.00    0.00   0.000
   0.000   0.000   0.000      0.      0.   0.000   0.000   0.000      0.      0.
"""

_GeoNet_file_format="""
visit: 
    http://info.geonet.org.nz/display/appdata/Accelerogram+Data+Filenames+and+Formats
    Blocks = {A, B, C, D}
    A: 
      Text header block
    B:
      int header block
    C:
      float header block
    D:
      float data blocks
"""

_event_origin="""\
year
month
day
hour
min
secx10
NZ_seism_obs_log_num
aftershock_num
buffer_start_time_year
buffer_start_time_month
""".split()

_source_info="""\
lat_deg
lat_min
lat_sec
lon_deg
lon_min
lon_sec
hypocentral_depth
centroid_depth
buffer_start_time_day
buffer_start_time_hour
""".split()

_site_info="""\
lat_deg
lat_min
lat_sec
lon_deg
lon_min
lon_sec
lon_axis_bearing
comp_dir
site_epicentral_bearing
epicentral_dist
""".split()

_sample_info="""\
num_digitised_samples
prepended_samples
appended_samples
acc_samples
vel_samples
disp_samples
unused1
unused2
buffer_start_time_minute
buffer_start_time_secx1000
""".split()

_line_21="""\
instrument_freq
ratio_critical_damping
film_speed
scale_factor_to_mm_on_film
inverse_sensitivity
deviation_angle_1
deviation_angle_2
factor_to_give_unit_of_mmperspers
timing_lamp_offset
time_of_common_timemark
""".split()

_line_22="""\
site_lat_deg
site_lon_deg
epicentral_lat_deg
epicentral_lon_deg
Ml
Ms
Mw
Mb
unused1
unused2
""".split()

_line_23="""\
digitised_duration
prepended_duration
appended_duration
record_duration
orig_sampling_interval
sampling_interval
unused1
unused2
unused3
local_g
""".split()

_line_24="""\
peak_unfiltered_accel
time_peak_unfiltered_accel
rms_unfiltered_accel
unused1
unused2
peak_filtered_accel
time_peak_filtered_accel
rms_filtered_accel
peak_horiz_accel
dominant_filtered_accel_freq
""".split()

_line_25="""\
peak_vel
time_peak_vel
rms_vel
peak_horiz_vel
dominant_vel_freq
peak_disp
time_peak_disp
rms_disp
peak_horiz_disp
dominant_disp_freq
""".split()

_line_26="""\
high_pass_cutoff_freq
high_pass_trans_freq
high_pass_rolloff_freq
unused1
unused2
low_pass_rolloff_freq
low_pass_trans_freq
low_pass_cutoff_freq
unused3
unused4
""".split()

import numpy as np
class FileComponent(object):
    
    def __init__(self):
        self.acc = None
        self.vel = None
        self.disp= None
        self.A_header = None
        self.lines = []
        self.angle = None
        self.delta_t = None
        self.time_delay = None
        self.B_header = {"event_origin":{}, 
                          "source_info":{},
                            "site_info":{},
                          "sample_info":{}}
        self.C_header = {"line_21":{},
                         "line_22":{},
                         "line_23":{},
                         "line_24":{},
                         "line_25":{},
                         "line_26":{}}
        #_GeoNet_file_format
    def extract(self,lines):
        self.A_header = lines[0:16]
        
        self.B_header["event_origin"].update(
        zip(_event_origin, np.asarray(lines[16].split(), dtype='i'))
        )
        self.B_header["source_info"].update(
        zip(_source_info, np.asarray(lines[17].split(), dtype='i'))
        )
        self.B_header["site_info"].update(
        zip(_site_info, np.asarray(lines[18].split(), dtype='i'))
        )
        self.B_header["sample_info"].update(
        zip(_sample_info, np.asarray(lines[19].split(), dtype='i'))
        )

        self.C_header["line_21"].update(
        zip(_line_21, np.asarray(lines[21-1].split(), dtype='f'))
        )
        self.C_header["line_22"].update(
        zip(_line_22, np.asarray(lines[22-1].split(), dtype='f'))
        )
        self.C_header["line_23"].update(
        zip(_line_23, np.asarray(lines[23-1].split(), dtype='f'))
        )
        self.C_header["line_24"].update(
        zip(_line_24, np.asarray(lines[24-1].split(), dtype='f'))
        )
        self.C_header["line_25"].update(
        zip(_line_25, np.asarray(lines[25-1].split(), dtype='f'))
        )
        self.C_header["line_26"].update(
        zip(_line_26, np.asarray(lines[26-1].split(), dtype='f'))
        )

        self.angle = float(self.B_header['site_info']['comp_dir'])
        self.time_delay = (self.B_header['sample_info']['buffer_start_time_minute']-
                           self.B_header['event_origin']['min'])*60000+\
                          (self.B_header['sample_info']['buffer_start_time_secx1000']-
                           self.B_header['event_origin']['secx10']*100)
        self.time_delay *=1e-3

        self.delta_t = float(self.C_header['line_23']['sampling_interval'])
        from math import ceil
        num_acc_lines = int(ceil(self.B_header["sample_info"]["acc_samples"]/10.))
        num_vel_lines = int(ceil(self.B_header["sample_info"]["vel_samples"]/10.))
        num_disp_lines= int(ceil(self.B_header["sample_info"]["disp_samples"]/10.))
        lines = lines[26:]
        self.acc = lines[0:num_acc_lines]
        lines    = lines[num_acc_lines:]
        self.vel = lines[0:num_vel_lines]
        lines    = lines[num_vel_lines:]
        self.disp= lines[0:num_disp_lines]
        lines    = lines[num_disp_lines:]

        self.lines.append(self.acc)
        if (len(self.vel) is not 0):
            self.lines.append(self.vel)
            self.lines.append(self.disp) 
        
        self.acc = [ np.asarray(_.split(),dtype='f') for _ in self.acc ]
        self.acc = np.concatenate(self.acc)

        if (len(self.vel) is not 0):
            self.vel = [ np.asarray(_.split(),dtype='f') for _ in self.vel ]
            self.vel = np.concatenate(self.vel)
            self.disp= [ np.asarray(_.split(),dtype='f') for _ in self.disp ]
            self.disp= np.concatenate(self.disp)



        return lines
    
    def __str__(self):
        for key, value in self.B_header.items():
            print("\n********** %s ***********\n" % key)
            print(value)
        
        print("\n******************************\n")
        return " ".join(self.A_header)

import os 
class GeoNet_File(object):
    
    def __init__(self, station_fileName, base_dir=os.getcwd(),vol=1):
        self.station_fileName = station_fileName
        self.vol = vol
        self.base_dir = base_dir
        self._parse()
    
    def _readFile(self):
        lines = []
        with open("/".join([self.base_dir, self.station_fileName]),'r') as f:
            lines = f.readlines()
        
        return lines
    
    def _parse(self):
        lines = self._readFile()
        self.comp_1st, self.comp_2nd, self.comp_up = (
                    FileComponent() for i in range(3))
        
        lines = self.comp_1st.extract(lines)
        lines = self.comp_2nd.extract(lines)
        lines = self.comp_up.extract(lines)
        
        assert (len(lines) is 0), "D'oh! Final list must be empty"
        return
    
if __name__ == "__main__":
    gf = GeoNet_File("/".join(["tests","data",
                     "20160214_001345_NBLC_20.V1A"]),vol=1)
    print(gf.comp_1st)

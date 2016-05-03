import numpy as np
import scipy.linalg
from scipy.integrate import cumtrapz
from butterworth import ButterWorth, butter_bandpass_filter

cos, sin, pi = (np.cos, np.sin, np.pi)


class SMD(object):
    """
    Strong Motion Data container
    """
    _template = dict.fromkeys(('BB','HF','LF'))
    
    def __init__(self, acc, dt=0.005, lowcut=0.1, highcut=50.,
                     fs=(1./0.005), ft=1., order=4, output=None): 

        self._acc = self._template.copy()
        self._vel = self._template.copy()
        self._disp= self._template.copy()
        
        self.acc = acc
        self.dt = dt
        
        self.lowcut=lowcut
        self.highcut=highcut
        self.fs = fs
        self.ft =ft
        self.order=order
        self.output=output
    
    @property
    def accBB(self):
        if self._acc['BB'] is None:
            self._acc['BB'] = butter_bandpass_filter(data=self.acc,
              lowcut=self.lowcut, highcut=self.highcut, fs=self.fs, 
                                                  order=self.order)
        return self._acc['BB']
        
    @property
    def velBB(self):
        if self._vel['BB'] is None:
            self._vel['BB'] = cumtrapz(y=self.accBB, dx=self.dt, initial=0.)

        return self._vel['BB']
        
    @property
    def dispBB(self):
        if self._disp['BB'] is None:
            self._disp['BB'] = cumtrapz(y=self.velBB, dx=self.dt, initial=0.)

        return self._disp['BB']

    @property
    def accHF(self):
        if self._acc['HF'] is None:
            self._acc['HF'] = ButterWorth.filter(data=self.accBB, 
                         btype='highpass', fs=self.fs, ft=self.ft, 
                             order=self.order, output=self.output)

        return self._acc['HF']

    @property
    def accLF(self):
        if self._acc['LF'] is None:
            self._acc['LF'] = ButterWorth.filter(data=self.accBB, 
                         btype='lowpass', fs=self.fs, ft=self.ft, 
                            order=self.order, output=self.output)

        return self._acc['LF']

    @property
    def velHF(self):
        if self._vel['HF'] is None:
            self._vel['HF'] = ButterWorth.filter(data=self.velBB, 
                        btype='highpass', fs=self.fs, ft=self.ft, 
                            order=self.order, output=self.output)

        return self._vel['HF']

    @property
    def velLF(self):
        if self._vel['LF'] is None:
            self._vel['LF'] = ButterWorth.filter(data=self.velBB, 
                         btype='lowpass', fs=self.fs, ft=self.ft, 
                            order=self.order, output=self.output)

        return self._vel['LF']

    @property
    def dispHF(self):
        if self._disp['HF'] is None:
            self._disp['HF'] = ButterWorth.filter(data=self.dispBB, 
                          btype='highpass', fs=self.fs, ft=self.ft, 
                              order=self.order, output=self.output)

        return self._disp['HF']

    @property
    def dispLF(self):
        if self._disp['LF'] is None:
            self._disp['LF'] = ButterWorth.filter(data=self.dispBB, 
                           btype='lowpass', fs=self.fs, ft=self.ft, 
                              order=self.order, output=self.output)

        return self._disp['LF']



class Process(object):

    def __init__(self, gf, automatic=True):
        
        self.gf = gf
        self.rot_angle = 0.
        self.delta_t = gf.comp_1st.delta_t

        self.acc_000 = None
        self.acc_090 = None
        self.acc_ver = gf.comp_up.acc.copy()

        assert gf.vol is 1, "Only vol1 data is processed."
        if automatic:
            self.rotate()

            SMD_kwargs = {'lowcut':0.1, 'highcut':50., 'fs':(1./0.005), 
                          'order':4, 'ft':1., 'output':None}

            SMD_kwargs['fs']=self.gf.comp_1st.C_header['line_21']['instrument_freq']
        
            self.comp_000 = SMD(self.acc_000, self.delta_t, **SMD_kwargs)
            self.comp_090 = SMD(self.acc_090, self.delta_t, **SMD_kwargs)
            self.comp_ver = SMD(self.acc_ver, self.delta_t, **SMD_kwargs)

    def rot_matrix(self,theta):
        """
            Rotating clockwise
            theta: rotation angle in degrees
        """
        theta = theta*pi/180.
        return np.array([[cos(theta), sin(theta)],
                        [-sin(theta), cos(theta)]])
    
    def rotate(self):
        """
                N
            W<-- -->E
                S
        
        GeoNet comp_1st and comp_2nd angles axis angles are measured
        from N. We rotate clock wise with angle theta where
            
            theta = 360 - comp_1st.angle
        
        [comp_000, comp_180] = Rot(theta) * [comp_1st, comp_2nd]
        Then perform a reflection in the y-axis
        comp_090 = -comp_180

        self.rot_angle: 
               theta = 360 - comp_1st.angle        
        acc_000:
                positive axis 1 parallel to N
        acc_090:
                positive axis 2 parallel to E
        
        """
        self.rot_angle = 360. - self.gf.comp_1st.angle
        R = self.rot_matrix(self.rot_angle)
        self.acc_000 = R[0,0] * self.gf.comp_1st.acc + R[0,1]*self.gf.comp_2nd.acc
        self.acc_090 = -(R[1,0] * self.gf.comp_1st.acc + R[1,1]*self.gf.comp_2nd.acc)
        return

    def save2disk(self,loc, stat_code):
        """
        The first two rows in .000, .090 and .ver files are
        row1: 
            Station code    comp   title(optional)
        e.g   CCCC           90     summed output (for broadband)
        row2:
            # data points delta(t) HR MN SEC  3 unused floats
        e.g	    22000      0.005   0. 0. -1.  0. 0. 0.  
        """
        print("\nSaving Rotated Broadband acceleration data for %s to: "
              % stat_code)
        print("%s \n" % loc)
        
        header_000 = stat_code + " 0 broadband\n"
        header_090 = stat_code + " 90 broadband\n"
        header_ver = stat_code + " ver broadband\n"
        time_delay = str(self.gf.comp_1st.time_delay)

        header_000 += " ".join(map(str, [self.comp_000.accBB.size,
                                         self.delta_t,
                                         "0. 0. "+time_delay+" 0. 0. 0."]
                                         ))
        header_090 += " ".join(map(str, [self.comp_090.accBB.size,
                                        self.delta_t,
                                        "0. 0. "+time_delay+" 0. 0. 0."]
                                        ))
        header_ver += " ".join(map(str, [self.comp_ver.accBB.size,
                                         self.delta_t,
                                         "0. 0. "+time_delay+" 0. 0. 0."]
                                         ))
        
        f = "/".join([loc,stat_code])
        np.savetxt(".".join([f,"000"]), self.comp_000.accBB,
                   header=header_000, comments="")
        np.savetxt(".".join([f,"090"]), self.comp_090.accBB,
                   header=header_090, comments="")
        np.savetxt(".".join([f,"ver"]), self.comp_ver.accBB,
                   header=header_ver, comments="")
        return


if __name__ == "__main__":
    from geoNet_file import GeoNet_File
    gf = GeoNet_File("/".join(["tests","data",
                     "20160214_001345_NBLC_20.V1A"]),vol=1)
    pgf = Process(gf)
    print(pgf.comp_000.accBB, pgf.comp_000.accLF, pgf.comp_000.accHF)

from geoNet_file import GeoNet_File
from process import Process
import butterworth as bw
import os


gf_wvz = GeoNet_File("20100903_163541_WVZ.V1A",
                     "/".join([os.getcwd(),"tests","2010Mw7pt2_data"]),
                     vol=1)
#normalize with g=9810. and get acc in cm/s^2
g = 9810. #mm/s^2
gf_wvz.comp_up.acc  *= g/gf_wvz.comp_up.C_header["line_23"]["local_g"]/10.
gf_wvz.comp_1st.acc *= -g/gf_wvz.comp_1st.C_header["line_23"]["local_g"]/10.
gf_wvz.comp_2nd.acc *= -g/gf_wvz.comp_2nd.C_header["line_23"]["local_g"]/10.

pgf_wvz = Process(gf_wvz)
#Read the same data as filtered by Brendon
floc = "/hpc/home/hnr12/ObservedGroundMotions/20100903163541/BB/"
fname_ver = "WVZ.ver"
import numpy as np
fver = open("/".join([floc,fname_ver]),'r')
vel_ver_BB = np.genfromtxt(fname=fver, dtype='float',
                           skip_header=2,skip_footer=0)
fname_000 = "WVZ.000"
import numpy as np
f000 = open("/".join([floc,fname_000]),'r')
vel_000_BB = np.genfromtxt(fname=f000, dtype='float',
                           skip_header=2,skip_footer=0)

fname_090 = "WVZ.090"
import numpy as np
f090 = open("/".join([floc,fname_090]),'r')
vel_090_BB = np.genfromtxt(fname=f090, dtype='float',
                           skip_header=2,skip_footer=0)


#last_line = f.readline()
#last_line = np.asarray([x for x in last_line.split()], dtype=float)
                        
fver.close()
f000.close()
f090.close()

#print("the last line was " %last_line)
vel_ver_BB = np.concatenate(vel_ver_BB)
vel_000_BB = np.concatenate(vel_000_BB)
vel_090_BB = np.concatenate(vel_090_BB)

vel_ver = pgf_wvz.comp_ver.velBB
vel_000 = pgf_wvz.comp_000.velBB
vel_090 = pgf_wvz.comp_090.velBB

from matplotlib import pylab as plt
delta_t = gf_wvz.comp_up.delta_t
t = np.arange(vel_ver_BB.size)*delta_t
#t = np.arange(vel_ver.size)*delta_t
fig = plt.figure()
fig.suptitle("Velocity WVZ.VER")
plt.subplot(121)
plt.plot(t,vel_ver_BB-vel_ver, 'k-',label='BB-Ah')
plt.legend(loc="best")

plt.subplot(122)
plt.plot(t, vel_ver_BB, 'k:',label='BB')
plt.plot(t, vel_ver, 'b--', label='Ah')
plt.legend(loc="best")

fig = plt.figure()
fig.suptitle("Velocity WVZ.000")
plt.subplot(121)
plt.plot(t,vel_000_BB-vel_000, 'k-',label='BB-Ah-000')
plt.legend(loc="best")

plt.subplot(122)
plt.plot(t, vel_000_BB, 'k:',label='BB-000')
plt.plot(t, vel_000, 'b--', label='Ah-000')
plt.legend(loc="best")

fig = plt.figure()
fig.suptitle("Velocity WVZ.090")
plt.subplot(121)
plt.plot(t,vel_090_BB-vel_090, 'k-',label='BB-Ah-090')
plt.legend(loc="best")

plt.subplot(122)
plt.plot(t, vel_090_BB, 'k:',label='BB-090')
plt.plot(t, vel_090, 'b--', label='Ah-090')
plt.legend(loc="best")


plt.show()


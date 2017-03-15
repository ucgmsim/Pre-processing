from geoNet_file import GeoNet_File
from process import Process
import butterworth as bw
import os
#from utilities import readGP
from utils import readGP
import numpy as np
gf_wvz = GeoNet_File("20100903_163541_WVZ.V1A",
                     "/".join([os.getcwd(),"tests","data"]),
                     vol=1)

pgf_wvz = Process(gf_wvz, lowcut=0.1)
#Read the same data as filtered by Brendon
floc = "/hpc/home/hnr12/ObservedGroundMotions/20100903163541/BB/"
vel_ver_BB , num_pts, dt, shift= readGP(floc, "WVZ.ver")
vel_000_BB , num_pts, dt, shift= readGP(floc, "WVZ.000")
vel_090_BB , num_pts, dt, shift= readGP(floc, "WVZ.090")

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


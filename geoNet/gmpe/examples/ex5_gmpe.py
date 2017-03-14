"""
Read stations_list.ll and find the Rrups for each SMS in it.
"""
import numpy as np
from itertools import izip
from geoNet.utils import read_statsll, get_empIM_v2
from geoNet.gmpe import readStationFile as rsf
from geoNet.gmpe.calculateGMPE import set_faultprop 



#Read the standard rupture file
srf_fname="/nesi/projects/nesi00213/RupModel/ahsan/2011Apr29_m4pt9/ptSource/Srf/2011Apr29_v2_m4pt9.srf"
FiniteFault = rsf.readSrfFile(srf_fname)
FiniteFault = rsf.Points_np(FiniteFault)

#set fault properties
faultprop = set_faultprop(Mw=4.9, rake=157., dip=80., Ztor=7.)


stats_floc="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017"
stats_fname="20110429_190804_eventStats_2017-01-10.ll"
event_stats=read_statsll(stats_floc, stats_fname)


periods=np.array([0.1, 0.5, 1.,  5., 10.])
Rrups = np.logspace(np.log10(1.),np.log10(250.),30)

empIM_values, empIM_sigmas =  get_empIM_v2([2,3],[1,2,3,4],faultprop)
print empIM_values
print "and"
print empIM_sigmas

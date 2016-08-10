# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 16:02:59 2015

@author: Richard Clare

parametersStation.py

Parameters to be used for post processing of for each station (calculation and plotting of the Intensity Measures and velocity time series).

Issues:
Should read runName from .e3d file not params.py? cleaner. used par_value from shared.py
Maybe should check to see if the IM .h5 file already exists
Commented section should be uncommented for Fitzroy, and lines 'to be removed' be removed
"""

# imports all the Classes used to set up the parms files. The default values are stored in this file
from parmStructure import *
import numpy as np
import os
import sys
sys.path.append(os.path.abspath(os.path.curdir))
ANALYSIS='/nesi/projects/nesi00213/groundMotionStationAnalysis'
import params as runparams

# flags control the logic in runPostProcess.py (what to calculate/plot)
# 0: none
# 1: observation only
# 2: simulation only
# 3: both
flag = Flag()
flag.calcIM = 3
flag.plotSeismogram = 3
flag.plotIM = 3
flag.plotTimeSeries = 3  # must be 0 or 3 only (neither or both)
flag.plotGMPE = 1        # 1: overlays GMPE on IM plot, 0 does not

#parameters
#GEN_ROCK_VS = 500, 
#KAPPA = 0.045,
#QFEXP = 0.6,
#SDROP = 50, 
#np.array([10, 13, 15. ,20, 25, 30, 33, 35, 37, 40, 45])


suffix = '_500_045_60_50' #'500_045_60_40'        
#
# where to find/place things
files = Files()
ObsDir_base = '/nesi/projects/nesi00213/ObservedGroundMotions/ahsan'
velLFObsDir = ObsDir_base + '/Mw4pt6_20101106_135203/Vol1/data/velLF'
velBBObsDir = ObsDir_base + '/Mw4pt6_20101106_135203/Vol1/data/velBB'

files.figsDir       = runparams.figures_dir
from params_base_bb import bb_veldir
files.velSimDir     = bb_veldir
files.velObsDir     = velBBObsDir

files.stationFile   = runparams.stat_file
files.srfFile       = runparams.srf_file
files.runName       = runparams.run_name#'2011Apr29_m4pt9'+suffix
del runparams

# siteClassFile used if plotting based on the site class of each station
files.siteClassFile = '/'.join([ANALYSIS,'Sept42010GmMetadata.ll'])
#replace all the '.' in the runName with 'p'
files.runName=files.runName.replace('.','p')    #no longer need this. To be removed
# IMs filename to save to, (in vel files dir)
files.outputFile = files.runName+'_IMs.h5'
# IMs code is able to calculate/save/plot
files.ImPossibilities = ['AI','CAV','Ds575','Ds595','PGA','PGV','pSA']
# additional params to save
files.additionals = ['name','DT','NT','period']
# component names
files.components = ['000','090','ver']
# extention of velocity component files
files.extensions = ['.000','.090','.ver']

# constants for the IM calculation, especially pSA
calcIM = CalcIM()
calcIM.IMsToCalculate = files.ImPossibilities #['AI','CAV','Ds575','Ds595','PGA','PGV','pSA']  #The IMs we actually calculate and save out (or plot)
calcIM.period = np.logspace(np.log10(0.01), np.log10(10), 100)  # used for the pSA calculation
calcIM.deltat = 0.005  # time step (s) for integration (interpolated if DT for the time series is larger)
calcIM.c = 0.05        # viscous damping ratio
calcIM.g = 981.        # gravity (cm s^-2)
calcIM.M = 1.          # fixed Mass values
calcIM.useCython=True  # Whether to use Cython for the pSA calculation. You need to compile rspectra first
calcIM.gamma = 0.5     # Chopra constant. Only used in Cython implementation. Not Python
calcIM.beta = 0.25     # Chopra constant. Only used in Cython implementation. Not Python

# runPlotVelocitySeismogramsWithDistance related
plotSeismogram = PlotSeismogram()
plotSeismogram.magnify = 3.
plotSeismogram.normalize = 1.

# runPlotObsSynIntensityMeasures related
plotIM = PlotIM()
plotIM.ImName = ['PGA', 'pSA']      # must be a list ([]), all are plotted.
# 0 = NS, 1 = EW, 2 = vertical, 3 = geometric mean
plotIM.component = [3]              # must be a list ([]), all are plotted.
plotIM.pSAPeriod = [0.1, 0.3, 0.5, 1.0, 3.0, 5., 10.0] # (s) must be a list ([]), all are plotted.
# whether to give each station a different symbol based on site Class
plotIM.siteClassMarkers = False
plotIM.unknown ='D'     #What class to plot if it is unknown or unspecified 
plotIM.pSAStations=[] #['CBGS','LINC','CMHS','HVSC','LPCC']   #The stations to plot the inidividual pSA vs period for. An empty list ([]) will plot everything
plotIM.saveData=True
# PlotObservedSimulatedTimeSeries related
plotTimeSeries = PlotTimeSeries()
# empty list ([]) will plot everything
plotTimeSeries.stationsToPlot = [] #['CBGS','LINC','CMHS','HVSC','LPCC']
plotTimeSeries.tmax = 100      # 100 useful default
plotTimeSeries.plotTimeSeriesAspectRatio = 1
plotTimeSeries.gapObsSim = 2   # 2 useful default (peaks at same time would touch)
plotTimeSeries.gapStations = 5 # 5 useful default (peaks at same time would be 1 unit apart)
plotTimeSeries.subdir = 'waveforms'
plotTimeSeries.nLimit=5
#Options for overlaying the GMPE in runPlotObsSynIntensityMeasures
plotGMPE = PlotGMPE()
plotGMPE.model = 1                          # 1 = Bradley2010 NZ-specific (PGA, SA) is only one currently coded
plotGMPE.DistMax = 100                      # maximum distance (km) for calculating the GMPE
plotGMPE.DistMin = 0.1                      # maximum distance (km) for calculating the GMPE
plotGMPE.nval = 51                          # no. of values used in generating the GMPE 
plotGMPE.SiteClassForPrediction = ['D']     # which site class (hence Vs30) to use for preidiciton. eg ['D', 'B']
plotGMPE.Vs30 = [1500, 800, 450, 250, 200]  # shear wave velocity (m/s) for different classes depending on SiteClassForPredcition. [A, B, C, D, E]
plotGMPE.Ztor = 6.#centroid depth                           # depth to top of coseismic rupture (km)
plotGMPE.Rtvz = 0                           # source-to-site distance in the Taupo volcanic zone (TVZ) (km)
plotGMPE.V30measured = 0                    # yes =1 (i.e. from Vs tests); no=0 (i.e. estimated from geology)
plotGMPE.Mw=4.6                             # moment tensor magnitude
plotGMPE.dip=87.                            # average dip angle (degrees)
plotGMPE.rake=157.                          # rake angle (degrees)

#Options for running a unit test
tests=Tests()
tests.IM_unit_test=False                            # Compare the IMs with a reference set (Darfield earthquake observations)
tests.Rrup_unit_test=False                          # Save out the rupture distances and compare with the reference values for this station file and srf file
tests.GMPE_unit_test=False                           # Compare the GMPEs for the Bradley2010 method with the Darfield earthquake observations
tests.IM_unit_test_filename=os.path.abspath(os.path.join(ANALYSIS,'Unit_Tests/IMs_unit_test.h5'))      # Name of the file with the IMs for unit test
tests.Rrup_unit_test_filename=os.path.abspath(os.path.join(ANALYSIS,'Unit_Tests/Rrups_unit_test.h5'))  # Name of the file with the Rrups for unit test
tests.GMPE_PGA_unit_test_filename=os.path.abspath(os.path.join(ANALYSIS,'Unit_Tests/GMPE_PGA_Bradley2010_unit_test.h5')) 
#Name of the file with the GMPEs for unit test for PGA
tests.GMPE_pSA_unit_test_filename=os.path.abspath(os.path.join(ANALYSIS,'Unit_Tests/GMPE_pSA1s_Bradley2010_unit_test.h5')) 
#Name of the file with the GMPEs for unit test for pSA at 1s

print files.__dict__.items()

"""
Created on Mon Mar 14 09:30:44 2016

@author: rmc84

Richard Clare

Purpose: calcualtes the GMPE values for a given IM & earthquake

Based on CompareGMPEs_.m (version 1.0 8 March 2010, Brendon Bradley)
All variable and function names have been retained wherever possible

All the redundant parts of the code to plot the GMPE on the IM plot have been removed.

This function is mostly just to put the variables from the parms into the class (matlab structures) siteprop and faultproperties.

Other GMPEs can be called from this function with the same format.

Issues:

Verification:
Matches matlab code to 9sig fig (IM_Rscaling, Rrupval, sigma_IM_Rscaling)

"""

import numpy as np
#from Bradley_2010_Sa import Bradley_2010_Sa
from geoNet.gmpe.Bradley_2010_Sa import Bradley_2010_Sa

class siteProperties:       #Class of site properties. initialize all attributes to None
    period=None             # '(-1),(0),(real variable)' period of vibration =-1->PGV; =0->PGA; >0->SA
    Rrup=None               # closest distance coseismic rupture (km)
    Rjb=None                # ???
    Rx=None                 #distance measured perpendicular to fault strike from surface projection of updip edge of the fault rupture (+ve in downdip dir) (km)
    Rtvz=None               # source-to-site distance in the Taupo volcanic zone (TVZ) (km)
    V30measured=None        # yes =1 (i.e. from Vs tests); no=0 (i.e. estimated from geology)
    V30=None                # shear wave velocity at 30m depth (m/s)
    Z1pt0=None              # depth to the 1.0km/s shear wave velocity horizon (optional, uses default relationship otherwise
    g=None                  # gravity (cm s^-2)

class faultProperties():    #Class of fault properties. initialize all attributes to None
    Mw=None                 # moment tensor magnitude
    rake=None               # rake angle (degrees)
    dip=None                # dip angle (degrees)
    Ztor=None               # depth to top of coseismic rupture (km)
                 

class set_faultprop(object):

    def __init__(self,Mw=None, rake=None, dip=None, Ztor=None):
        """
        In one line set the fault parameters
        """

        #self.faultprop=faultProperties()
        #self.faultprop.Mw=Mw
        #self.faultprop.rake=rake
        #self.faultprop.dip=dip
        #self.faultprop.Ztor=Ztor 

        self.Mw=Mw
        self.rake=rake
        self.dip=dip
        self.Ztor=Ztor 

        return

class set_siteprop(object):

    def __init__(self,faultprop,
        period=None, Rrup=None, Rjb=None, Rx=None, Rtvz=None, V30measured=None,
        V30=None, Z1pt0=None, g=981):
        """
        g in cm/s^2
        Period = -1 for PGV, 0 for PGA, actual numerical value for PSA
        Default choice in Richard's code used is V30=250. #site 'D' classification
        """

#        self.siteprop=siteProperties()
#        self.siteprop.period=period
#        self.siteprop.g=g      
#        self.siteprop.Rtvz=Rtvz     
#        self.siteprop.V30measured=V30measured 
#        #self.siteprop.V30=250. #site 'D' classification
#        self.siteprop.V30=V30 
#        self.siteprop.Rrup=Rrup
#        #Definition used below by Richard does not match Rjb definition
#        self.siteprop.Rjb=np.sqrt(np.max((0,Rrup**2-faultprop.Ztor**2)))
#        #We leave Rjb unchanged to match Ricard's code
#        #self.siteprop.Rjb=Rjb
#        self.siteprop.Rx=-self.siteprop.Rjb
#
        self.period=period
        self.Z1pt0=Z1pt0
        self.g=g      
        self.Rtvz=Rtvz     
        self.V30measured=V30measured 
        #self.V30=250. #site 'D' classification
        self.V30=V30 
        self.Rrup=Rrup
        #Definition used below by Richard does not match Rjb definition
        self.Rjb=np.sqrt(np.max((0,Rrup**2-faultprop.Ztor**2)))
        #We leave Rjb unchanged to match Ricard's code
        #self.Rjb=Rjb
        self.Rx=-self.Rjb
        return 


def calculateGMPE(parms,ImName,period):     ##need to pass the ImName and the period as we loop over these

    #Do some error checking to make sure that the model is valid

    if parms.plotGMPE.model!=1:
        print ('Error: The selected GMPE (model %d) is not supported. Only Model=1 Bradley_2010 is supported.') %(parms.plotGMPE.model)
        print ('Not calculating the GMPE. Returning from function')
        raise ValueError
        return

    #For the model, check that the IM is supported
    if not((ImName=='PGA')|(ImName=='PGV')|(ImName=='pSA')):
        print ('Error: The selected GMPE (model %d) does not support IM name %s') %(parms.plotGMPE.model,ImName)
        print ('Not calculating the GMPE for %s. Returning from function') %ImName
        raise ValueError
        return

    #set up the class faultprop. Same structure as the matlab code
    faultprop=faultProperties()
    faultprop.Mw=parms.plotGMPE.Mw
    faultprop.rake=parms.plotGMPE.rake
    faultprop.dip=parms.plotGMPE.dip
    faultprop.Ztor=parms.plotGMPE.Ztor 

    #set up the class siteprop. Same structure as the matlab code
    siteprop=siteProperties()
    if ImName=='PGA':
        siteprop.period=0
    elif ImName=='PGV':
        siteprop.period=-1
    elif ImName=='pSA':
        siteprop.period=period

    siteprop.g=parms.calcIM.g      
    siteprop.Rtvz=parms.plotGMPE.Rtvz     
    siteprop.V30measured=parms.plotGMPE.V30measured 

    #Rrup values used for the GMPE calculation
    Rrupval=np.exp(np.linspace(np.log(parms.plotGMPE.DistMin),np.log(parms.plotGMPE.DistMax),parms.plotGMPE.nVal))     

    Vs30=[]     #The list of Vs30 values we are going to calculate the GMPE for
    for i in range(len(parms.plotGMPE.SiteClassForPrediction)):
        if parms.plotGMPE.SiteClassForPrediction[i]=='A':
            Vs30.append(parms.plotGMPE.Vs30[0])
        elif parms.plotGMPE.SiteClassForPrediction[i]=='B':
            Vs30.append(parms.plotGMPE.Vs30[1])
        elif parms.plotGMPE.SiteClassForPrediction[i]=='C':
            Vs30.append(parms.plotGMPE.Vs30[2])
        elif parms.plotGMPE.SiteClassForPrediction[i]=='D':
            Vs30.append(parms.plotGMPE.Vs30[3])
        elif parms.plotGMPE.SiteClassForPrediction[i]=='E':
            Vs30.append(parms.plotGMPE.Vs30[4])
        else:
            print (parms.plotGMPE.SiteClassForPrediction[i], 'is an invalid siteClass for Prediciton. Must be A, B, C, D, E.')

    #initialise the outputs as arrays of zeros
    IM_Rscaling=np.zeros((len(Vs30),parms.plotGMPE.nVal))
    sigma_IM_Rscaling=np.zeros((len(Vs30),parms.plotGMPE.nVal,3))
    
    for j in range(len(Vs30)):  #loop over the different Vs30 values chosen
        siteprop.V30=Vs30[j] 
        for i in range(len(Rrupval)):   #loop over the different rupture distances
            siteprop.Rrup=Rrupval[i]
            siteprop.Rjb=np.sqrt(np.max((0,Rrupval[i]**2-faultprop.Ztor**2)))
            siteprop.Rx=-siteprop.Rjb            
            if parms.plotGMPE.model==1:     #Bradley2010 NZ-specific. No error checking here as we check at start of file                
                IM_Rscaling[j,i], sigma_IM_Rscaling[j,i,:]=Bradley_2010_Sa(siteprop,faultprop)                    

    return Rrupval, Vs30, IM_Rscaling, sigma_IM_Rscaling


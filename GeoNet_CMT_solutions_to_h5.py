"""
GeoNet_CMT_solutions_to_h5.py
Richard Clare
Created 03/05/2016

Script to 

1) read in GeoNet_CMT_solutions.csv and save as .h5 file.

2) load .h5 file and sort depending on magnitude and lat and lon

header is:
EVENT_ID,Date,Latitude,Longitude,str,dp,rake,str,dp,rake,ML,Mw,Mo,CD,NS,DC,Mxx,Mxy,Mxz,Myy,Myz,Mzz,VR,Tva,Tpl,Taz,Nva,Npl,Naz,Pva,Ppl,Paz


"""

import h5rw
import sys
import string
import numpy as np

filename='GeoNet_CMT_solutions.csv'

out_filename='GeoNet_CMT_solutions.h5'

try:
    fid=open(filename)
    my_data = np.genfromtxt(filename, delimiter=',', skip_footer=7, skip_header=1)
except IOError:
    print 'Could not open' + filename
    sys.exit(0)

temp=fid.readlines() #reads line at a time - necessary for header
header=temp[0]
fields=header.split(',')
header.split(',')[-1].split('\r')[0] #remove the end of line from the last one

#Put the data in a dictionary
quakes={}
for i in range(len(fields)):
    quakes[fields[i]]=my_data[:,i]

#save out the unsorted dictionary as .h5
h5rw.h5write(out_filename,quakes)

#############################################################################
#2) choose the earthquakes we want based on Mw, lat and lon

#set the wanted minimum and maximum latitudes
#these are for South Island only
wanted_lat=[-46.7,-40.4]    #min and max respectively
wanted_lon=[166.0,174.5]    #min and max respectively

min_Mw=5.0   #choose all quakes above this value

chosen=[]   #indices of the earthquakes we are choosing

#loop over all the quakes
for j in range(len(quakes[fields[0]])):
    lat=quakes['Latitude'][j]
    lon=quakes['Longitude'][j]
    Mw=quakes['Mw'][j]
    if (lat>wanted_lat[0])and(lat<wanted_lat[1])and(lon>wanted_lon[0])and(lon<wanted_lon[1])and(Mw>min_Mw):
        chosen.append(j)

#now create a dictionary with just the chosen quakes
chosen_quakes={}
for i in range(len(fields)):
    chosen_quakes[fields[i]]=quakes[fields[i]][chosen]

chosen_out_filename='chosen_'+out_filename

h5rw.h5write(chosen_out_filename,chosen_quakes)



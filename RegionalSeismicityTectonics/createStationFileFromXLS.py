# -*- coding: utf-8 -*-
"""
createStationFileFromXLS.py

Created on Thu Mar 17 10:12:03 2016

Script to open the GeoNetStations.xls spreadsheet and choose the stations by latitude and longitude to a stations file .ll

.ll file Tested on the Python Post-processing and gave the same output.
Not tested on the Emod3D code.

Issues:
1)I have the xlrd module on my local machine but not on Beatrice
2)Assumes a square area from lat,lon at the moment

@author: rmc84
"""

import numpy as np
import xlrd

xlsName='GeoNetStationCoordinates.xlsx'
stationFilename='southislandstations_v1.ll'#'mycantstations.ll'

wanted=['Station_Code', 'WGS84_LAT','WGS84_LON']
wanted_index=[]

wanted_lat=[-46.7,-40.4]    #min and max respectively
wanted_lon=[166.0,174.5]    #min and max respectively

#1) Open the .xls file 
xl_workbook = xlrd.open_workbook(xlsName)
sheet_names = xl_workbook.sheet_names()
xl_sheet = xl_workbook.sheet_by_index(0)
num_cols = xl_sheet.ncols   # Number of columns
num_rows = xl_sheet.nrows   # Number of rows

#2) Define the rows and columns that we want to save to the stations file
row = xl_sheet.row(0)  # 1st row
#loop over first row to find the wanted columns
for idx, cell_obj in enumerate(row):
    for idx2 in range(len(wanted)):
        if cell_obj.value==wanted[idx2]:
            wanted_index.append(idx)

#initialise the lists
chosen=[]       #list of indices of the rows we want to extract
lats=[]
lons=[]
codes=[]

""" Old code based on Station Code SC
chosen_network_codes='SC'
col = xl_sheet.col(0)   #This give the Column of Network Codes - colum 0
#loop over the columns to select the stations
for idx, cell_obj in enumerate(col):
    #get the rows of the chosen stations
    chosen.append(cell_obj.value in chosen_network_codes)

#loop over the rows and get the code, lat,lon from the chosen rows

for i in np.arange(num_rows):
    if chosen[i]==True:
        row = xl_sheet.row(i)
        #loop over the row elements
        for idx, cell_obj in enumerate(row):
            if idx==wanted_index[0]:            
                codes.append(cell_obj.value)
            if idx==wanted_index[1]:
                lats.append(cell_obj.value)
            if idx==wanted_index[2]:
                lons.append(cell_obj.value)
"""

#loop over to see if the lat and lon lie inside the box
for i in np.arange(num_rows):
    row = xl_sheet.row(i)
    #loop over the row elements
    for idx, cell_obj in enumerate(row):
        if idx==wanted_index[1]:
            lat=cell_obj.value
        if idx==wanted_index[2]:
            lon=cell_obj.value
    if (lat>wanted_lat[0])&(lat<wanted_lat[1])&(lon>wanted_lon[0])&(lon<wanted_lon[1]):
        chosen.append(i)

#loop over the rows and get the code, lat,lon from the chosen rows        
for i in np.arange(len(chosen)):
    row = xl_sheet.row(chosen[i])
    #loop over the row elements
    for idx, cell_obj in enumerate(row):
       if idx==wanted_index[0]:            
            codes.append(cell_obj.value)
        if idx==wanted_index[1]:
            lats.append(cell_obj.value)
        if idx==wanted_index[2]:
            lons.append(cell_obj.value)

#3) Write the required data to a binary file
fid=open(stationFilename,'w')
for i in np.arange(len(codes)):
    line='%.4f ' %lons[i] + '%.4f ' %lats[i] + codes[i] +'\n'
    fid.write(line)
    
fid.close()                

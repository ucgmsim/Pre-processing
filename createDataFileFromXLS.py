# -*- coding: utf-8 -*-
"""
createDataFileFromXLS.py

Created on Thu Mar 17 10:12:03 2016

Script to open the metadata .xls spreadsheet and choose the stations and data required
and write these to a txt file

We want:
StationCode
Lon
Lat
SiteClass
Vs30
VsMeasured
z1.0
Rjb
Rrup

@author: rmc84
"""

import numpy as np
import xlrd

xlsName='Sept42010GmMetadata.xlsx'
outputFilename=xlsName.split('.')[0]+'.ll'

'southislandstations_v1.ll'#'mycantstations.ll'

IMwanted='PGA(g)'

wanted=['Site code', 'Site class', 'Site Lat','Site Lon','Vs30measured','Vs30 (m/s)','Z1.0 (m)','Dept to Top Of Fault Rupture Model','Rjb (km)', 'Rrup (km)',IMwanted]
wanted_index=[]

wanted_lat=[-46.7,-40.4]    #min and max respectively
wanted_lon=[166.0,174.5]    #min and max respectively

#1) Open the .xls file 
xl_workbook = xlrd.open_workbook(xlsName)
sheet_names = xl_workbook.sheet_names()
xl_sheet = xl_workbook.sheet_by_index(0)
num_cols = xl_sheet.ncols   # Number of columns
num_rows = xl_sheet.nrows   # Number of columns

#2) Define the rows and columns that we want to save to the stations file
row = xl_sheet.row(0)  # 1st row
#loop over first row to find 
for idx, cell_obj in enumerate(row):
    for idx2 in range(len(wanted)):
        if cell_obj.value==wanted[idx2]:
            wanted_index.append(idx)

#2D array (list of lists) as both numbers and strings
all_data=[]

#loop over the rows and get all the wanted data        
for i in np.arange(1,num_rows):     #ignore the first row
    row = xl_sheet.row(i)
    all_data.append([])
    #loop over the row elements
    for j, cell_obj in enumerate(row):    
        for k, wi  in enumerate(wanted_index):        
            if j==wi:
                all_data[i-1].append(cell_obj.value)
        
#3) Write the required data to a binary file

fid=open(outputFilename,'w')
for i in np.arange(len(all_data)):
    line=''#initialize empty string
    for j in np.arange(len(all_data[0])):
        if type(all_data[i][j]) is float:
            line=line+ '%.6f ' %all_data[i][j]
        else:   #assume it's a string
            line=line+all_data[i][j]+' '
    line=line+'\n'            
    fid.write(line)

fid.close()        

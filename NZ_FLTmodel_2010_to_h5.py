"""
NZ_FLTmodel_20110_to_h5.py
Richard Clare
Created 03/05/2016

Script to 

1) read in NZ_FLTmodel_2010_Test.txt and save as .h5 file.

2) load .h5 file and sort dep

Need to remove the \n from each line

Row 1: FaultName 
Row 2: TectonicType , FaultType 
Row 3: LengthMean , LengthSigma (km) 
Row 4: DipMean , DipSigma (deg) 
Row 5: DipDir 
Row 6: Rake (deg) 
Row 7: RupDepthMean , RupDepthSigma (km) 
Row 8: RupTopMean, RupTopMin RupTopMax  (km) 
Row 9: SlipRateMean , SlipRateSigma (mm/yr) 
Row 10: CouplingCoeff , CouplingCoeffSigma (mm/yr) 
Row 11: MwMedian , RecurIntMedian  (yr) 
Row 12: Num Locations on Fault Surface 
Row 13+: Location Coordinates (Long, Lat) 


"""

import h5rw
import sys
import string
import numpy as np

filename='NZ_FLTmodel_2010_Test.txt'

out_filename='NZ_FLTmodel_2010_Test.h5'

try:
    fid=open(filename)
except IOError:
    print 'Could not open' + filename
    sys.exit(0)

temp=fid.readlines() #reads line at a time - necessary for headertemp=fid.readlines() #reads line at a time - necessary for header

#should store as a dictionary of lists to save as .h5
nlines=len(temp)
nlines_header=14

fields=['FaultName', 'TectonicType', 'FaultType', 'LengthMean', 'LengthSigma', 'DipMean', 'DipSigma', 'DipDir', 'Rake', 'RupDepthMean', 'RupDepthSigma', \
'RupTopMean', 'RupTopMin', 'RupTopMax', 'SlipRateMean', 'SlipRateSigma', 'CouplingCoeff', 'CouplingCoeffSigma', 'MwMedian', 'RecurIntMedian', 'NumLocations', \
'LocationCoordinates']

nFields=len(fields)

faults={}

#loop over all the fields initializing the list
for i in range(len(fields)):
    faults[fields[i]]=[]

#lines 0 to 14 are header - ignore

#loop over all the faults and all the fields
fault_counter=0
line_counter=nlines_header
while line_counter < nlines:    

    faults[fields[0]].append(temp[line_counter+1].split(' \n')[0])

    faults[fields[1]].append(temp[line_counter+2].split(' ')[0])
    faults[fields[2]].append(temp[line_counter+2].split(' ')[1])

    faults[fields[3]].append(np.float(temp[line_counter+3].strip().split(' ')[0]))
    faults[fields[4]].append(np.float(temp[line_counter+3].strip().split(' ')[-1]))

    faults[fields[5]].append(np.float(temp[line_counter+4].strip().split(' ')[0]))
    faults[fields[6]].append(np.float(temp[line_counter+4].strip().split(' ')[-1]))

    faults[fields[7]].append(np.float(temp[line_counter+5].strip().split(' \n')[0]))

    faults[fields[8]].append(np.float(temp[line_counter+6].strip().split(' \n')[0]))

    faults[fields[9]].append(np.float(temp[line_counter+7].strip().split(' ')[0]))
    faults[fields[10]].append(np.float(temp[line_counter+7].strip().split(' ')[-1]))

    faults[fields[11]].append(np.float(temp[line_counter+8].strip().split('   ')[0]))
    faults[fields[12]].append(np.float(temp[line_counter+8].strip().split('   ')[1].strip()))       
    faults[fields[13]].append(np.float(temp[line_counter+8].strip().split('   ')[2].strip()))

    faults[fields[14]].append(np.float(temp[line_counter+9].strip().split(' ')[0]))
    faults[fields[15]].append(np.float(temp[line_counter+9].strip().split(' ')[-1]))

    faults[fields[16]].append(np.float(temp[line_counter+10].strip().split(' ')[0]))
    faults[fields[17]].append(np.float(temp[line_counter+10].strip().split(' ')[-1]))

    faults[fields[18]].append(np.float(temp[line_counter+11].strip().split(' ')[0]))
    faults[fields[19]].append(np.float(temp[line_counter+11].strip().split(' ')[-1]))

    faults[fields[20]].append(np.int(temp[line_counter+12].strip().split(' \n')[0]))

    coords=[]
    for j in range(faults['NumLocations'][fault_counter]):
        lat=np.float(temp[line_counter+13+j].strip().split(' ')[0])
        lon=np.float(temp[line_counter+13+j].strip().split(' ')[-1])
        location=[lat,lon]
        coords.append(location)
    
    faults[fields[21]].append(coords)  
    
    line_counter=line_counter+13+faults['NumLocations'][fault_counter]
    fault_counter=fault_counter+1

#save out the dict to .h5
h5rw.h5write(out_filename,faults)

#############################################################################################
#Read in .h5 and sort based on lat and lon
#Different assumptions possible: all or part of the fault in the lat lon area

faults=h5rw.h5read(out_filename)

#set the wanted minimum and maximum latitudes
#these are for South Island only
wanted_lat=[-46.7,-40.4]    #min and max respectively
wanted_lon=[166.0,174.5]    #min and max respectively

#loop over all faults and check if it is wanted
nFaults=len(faults['FaultName'])

#list of indices of the faults in lat,lon range
chosen=[]

for m in range(nFaults):
    #also need to loop over all of NumLocations for each fault    
    for n in range(faults['NumLocations'][m]):        
        lat=faults['LocationCoordinates'][m][n][1]
        lon=faults['LocationCoordinates'][m][n][0]
        if (lat>wanted_lat[0])and(lat<wanted_lat[1])and(lon>wanted_lon[0])and(lon<wanted_lon[1]):
            chosen.append(m)
            break  #exit for loop if one of the locations is true


#create a new dictionary with just the chosen faults
chosen_faults={}

#loop over all the fields initializing the list
for i in range(len(fields)):
    chosen_faults[fields[i]]=[]

#loop over the fields
for ii in range(nFields):
    for jj in range(len(chosen)):
        chosen_faults[fields[ii]].append(faults[fields[ii]][chosen[jj]])

chosen_out_filename='chosen_'+out_filename

#now save out the chosen faults as .h5
h5rw.h5write(chosen_out_filename,chosen_faults)

#save out the fault as an .h5 file

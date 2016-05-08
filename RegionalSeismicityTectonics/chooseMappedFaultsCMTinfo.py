"""
chooseMappedFaultsCMTinfo.py
Richard Clare
Created 04/05/2016

Script to:

1) read in NZ_FLTmodel_2010_Test.h5 (mapped faults) and GeoNet_CMT_solutions.h5 (CMT info)

2) choose the mapped faults and CMT info based on lat, lon and Mw

3) save the chosen mapped faults and CMT info to .h5 again

"""

import h5rw
import sys
import string
import numpy as np

def NZ_FLTmodel_2010_to_h5():
    """Function to save the file NZ_FLTmodel_2010_Test.txt to a .h5 for easier manipulation.
        Only needs to be run whhen NZ_FLTmodel_2010_Test.txt is updated.

        No inputs. No outputs.

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

    mappedFaults={}

    #loop over all the fields initializing the list
    for i in range(len(fields)):
        mappedFaults[fields[i]]=[]

    #lines 0 to 14 are header - ignore

    #loop over all the faults and all the fields
    fault_counter=0
    line_counter=nlines_header
    while line_counter < nlines:    

        mappedFaults[fields[0]].append(temp[line_counter+1].split(' \n')[0])

        mappedFaults[fields[1]].append(temp[line_counter+2].split(' ')[0])
        mappedFaults[fields[2]].append(temp[line_counter+2].split(' ')[1])

        mappedFaults[fields[3]].append(np.float(temp[line_counter+3].strip().split(' ')[0]))
        mappedFaults[fields[4]].append(np.float(temp[line_counter+3].strip().split(' ')[-1]))

        mappedFaults[fields[5]].append(np.float(temp[line_counter+4].strip().split(' ')[0]))
        mappedFaults[fields[6]].append(np.float(temp[line_counter+4].strip().split(' ')[-1]))

        mappedFaults[fields[7]].append(np.float(temp[line_counter+5].strip().split(' \n')[0]))

        mappedFaults[fields[8]].append(np.float(temp[line_counter+6].strip().split(' \n')[0]))

        mappedFaults[fields[9]].append(np.float(temp[line_counter+7].strip().split(' ')[0]))
        mappedFaults[fields[10]].append(np.float(temp[line_counter+7].strip().split(' ')[-1]))

        mappedFaults[fields[11]].append(np.float(temp[line_counter+8].strip().split('   ')[0]))
        mappedFaults[fields[12]].append(np.float(temp[line_counter+8].strip().split('   ')[1].strip()))       
        mappedFaults[fields[13]].append(np.float(temp[line_counter+8].strip().split('   ')[2].strip()))

        mappedFaults[fields[14]].append(np.float(temp[line_counter+9].strip().split(' ')[0]))
        mappedFaults[fields[15]].append(np.float(temp[line_counter+9].strip().split(' ')[-1]))

        mappedFaults[fields[16]].append(np.float(temp[line_counter+10].strip().split(' ')[0]))
        mappedFaults[fields[17]].append(np.float(temp[line_counter+10].strip().split(' ')[-1]))

        mappedFaults[fields[18]].append(np.float(temp[line_counter+11].strip().split(' ')[0]))
        mappedFaults[fields[19]].append(np.float(temp[line_counter+11].strip().split(' ')[-1]))

        mappedFaults[fields[20]].append(np.int(temp[line_counter+12].strip().split(' \n')[0]))

        coords=[]
        for j in range(mappedFaults['NumLocations'][fault_counter]):
            lat=np.float(temp[line_counter+13+j].strip().split(' ')[0])
            lon=np.float(temp[line_counter+13+j].strip().split(' ')[-1])
            location=[lat,lon]
            coords.append(location)
        
        mappedFaults[fields[21]].append(coords)  
        
        line_counter=line_counter+13+mappedFaults['NumLocations'][fault_counter]
        fault_counter=fault_counter+1

    #save out the dict to .h5
    h5rw.h5write(out_filename,mappedFaults)


def GeoNet_CMT_solutions_to_h5():
    """

    Read in GeoNet_CMT_solutions.csv and save as .h5 file for easier maniuplation.
    Only needs to be run when GeoNet_CMT_solutions.csv is updated

    No inputs. No output.

    header is:
    EVENT_ID,Date,Latitude,Longitude,str,dp,rake,str,dp,rake,ML,Mw,Mo,CD,NS,DC,Mxx,Mxy,Mxz,Myy,Myz,Mzz,VR,Tva,Tpl,Taz,Nva,Npl,Naz,Pva,Ppl,Paz

    """

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
    CMTinfo={}
    for i in range(len(fields)):
        CMTinfo[fields[i]]=my_data[:,i]

    #save out the unsorted dictionary as .h5
    h5rw.h5write(out_filename,CMTinfo)


def chooseFaultsQuakes(min_Mw,wanted_lat,wanted_lon,faults_filename,quakes_filename,chosen_faults_filename,chosen_quakes_filename):
    """Function to read in the .h5 files of all faults and all earthquakes, and choose those in a specified latitude,longitude range and above a certain Mw threshold.
        Saves the chosen quakes and faults to .h5 files in the same dir.
    Inputs: wanted_lat -            list of [minimum, maximum] latitude for chosen area
            wanted_lon -            list of [minimum, maximum] longitude for chosen area
            min_Mw      -           minimum Mw (magnitude) threshold of earthquakes to choose
            mappedFaults_filename -        .h5 filename of all mapped faults
            CMTinfo_filename -       .h5 filename of all CMT info
            chosen_mappedFaults_filename - .h5 filename of the chosen mapped faults
            chosen_CMTinfo_filename - .h5 filename of the chosen CMT info

    Outputs: None.
        """

    #############################################################################################
    #1)Read in faults .h5 and sort based on lat and lon
    #Different assumptions possible: all or part of the fault in the lat lon area

    mappedFaults=h5rw.h5read(mappedFaults_filename)

    fields=mappedFaults.keys()
    nFields=len(fields)

    #loop over all faults and check if it is wanted
    nFaults=len(mappedFaults['FaultName'])

    #list of indices of the faults in lat,lon range
    chosen_mappedFault_indices=[]

    for m in range(nFaults):
        #also need to loop over all of NumLocations for each fault    
        for n in range(mappedFaults['NumLocations'][m]):        
            lat=mappedFaults['LocationCoordinates'][m][n][1]
            lon=mappedFaults['LocationCoordinates'][m][n][0]
            if (lat>wanted_lat[0])and(lat<wanted_lat[1])and(lon>wanted_lon[0])and(lon<wanted_lon[1]):
                chosen_mappedFault_indices.append(m)
                break  #exit for loop if one of the locations is true

    #create a new dictionary with just the chosen faults
    chosen_mappedFaults={}

    #loop over all the fields initializing the list
    for i in range(len(fields)):
        chosen_mappedFaults[fields[i]]=[]

    #loop over the fields
    for ii in range(nFields):
        for jj in range(len(chosen_mappedFault_indices)):
            chosen_mappedFaults[fields[ii]].append(mappedFaults[fields[ii]][chosen_mappedFault_indices[jj]])
   
    #now save out the chosen faults as .h5
    h5rw.h5write(chosen_mappedFaults_filename,chosen_mappedFaults)

    #############################################################################################
    #2)Read in quakes .h5 and sort based on lat and lon
    
    CMTinfo=h5rw.h5read(CMTinfo_filename)

    fields=CMTinfo.keys()

    chosen_CMTinfo_indices=[]   #indices of the earthquakes we are choosing

    #loop over all the quakes
    for j in range(len(CMTinfo['Mw'])):
        lat=CMTinfo['Latitude'][j]
        lon=CMTinfo['Longitude'][j]
        Mw=CMTinfo['Mw'][j]
        if (lat>wanted_lat[0])and(lat<wanted_lat[1])and(lon>wanted_lon[0])and(lon<wanted_lon[1])and(Mw>min_Mw):
            chosen_CMTinfo_indices.append(j)

    #now create a dictionary with just the chosen quakes
    chosen_CMTinfo={}
    for i in range(len(fields)):
        chosen_CMTinfo[fields[i]]=CMTinfo[fields[i]][chosen_CMTinfo_indices]
    
    h5rw.h5write(chosen_CMTinfo_filename,chosen_CMTinfo)


###########################################################################################################
#main

#only run these 2 lines if the CMTinfo and mappedFault files have been updated
NZ_FLTmodel_2010_to_h5()
GeoNet_CMT_solutions_to_h5()

#set the wanted minimum and maximum latitudes
#these are for South Island only
lat=[-46.7,-40.4]    #min and max respectively
lon=[166.0,174.5]    #min and max respectively

Mw=5.0      #set the minmum Magnitude threshold

mappedFaults_filename='NZ_FLTmodel_2010_Test.h5'    #name of the .h5 file containing all mappedFaults
CMTinfo_filename='GeoNet_CMT_solutions.h5'          #name of the .h5 file containing all CMTinof

chosen_mappedFaults_filename='chosen_NZ_FLTmodel_2010_Test.h5'      #name of the .h5 file where we save the chosen mappedFaults
chosen_CMTinfo_filename='chosen_GeoNet_CMT_solutions.h5'            #name of the .h5 file where we save the chosen CMTinfo

chooseFaultsQuakes(Mw,lat,lon,mappedFaults_filename,CMTinfo_filename,chosen_mappedFaults_filename,chosen_CMTinfo_filename)

#to load in the chosen faults and quakes
chosen_mappedFaults=h5rw.h5read(chosen_mappedFaults_filename)
chosen_CMTinfo=h5rw.h5read(chosen_CMTinfo_filename)




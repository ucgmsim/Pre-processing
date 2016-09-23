import os.path
import subprocess
binpath="/nesi/projects/nesi00213/Pre-processing/bin_process/devel"

def execute_cmd(cmd):
    p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outputList = p.stdout.readlines()
    output = ''.join(outputList)
    return output


if __name__=='__main__':

    ll=open('/nesi/projects/nesi00213/StationInfo/fd_sinz01-h0.400.ll')
    statcords=open('/nesi/projects/nesi00213/StationInfo/fd_sinz01-h0.400.statcords')
    statcords_lines=statcords.readlines()
    ll_lines=ll.readlines()
    assert(len(ll_lines)+1 == len(statcords_lines))

    f_processed=open('processed_stations.txt')
    l_processed=f_processed.readlines()
    l_processed=[x.strip('\n') for x in l_processed]

    station_collection={}

    for i in range(len(ll_lines)):
        columns=[x for x in ll_lines[i].strip('\n').split(' ') if x!='']
        station_collection[columns[2]]={'ll':(columns[0],columns[1]),'xy':None}


    for i in range(1,len(statcords_lines)):
        columns=[x for x in statcords_lines[i].strip('\n').split(' ') if x!='']
        station_collection[columns[3]]['xy']=(columns[0],columns[1])


    for actual_station in station_collection:
        if actual_station in l_processed:
            continue

        ll=station_collection[actual_station]['ll']
        xy=station_collection[actual_station]['xy']

        res= execute_cmd("python %s/hdf2ascii.py %s %s"%(binpath,ll[0],ll[1]))
        print "Station: %s  %s %s" %(actual_station,ll[0],ll[1])
        print "Expected XY: %s %s" %(xy[0],xy[1])
        print res
        zipfile=res.split('\n')[-2] #zipfile name
        virtual_station =zipfile.split('.zip')[0]
        res=execute_cmd("python %s/statgrid_validation.py %s %s"%(binpath,virtual_station,actual_station))
        print res



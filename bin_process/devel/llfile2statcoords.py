import tools
from params import *

filename='cantstations'
llfile=open(filename+'.ll','r')
statcoordsfile =open(filename+'.statcords','w')

ll_lines=llfile.readlines()
ll_inputs = []
for l in ll_lines:
   parts = l.strip('\n').split(' ')
   parts = [x for x in parts if x!='']
   ll_inputs.append((parts[2],parts[0],parts[1]))

statcoordsfile.write('%s\n'%len(ll_inputs))

for name, lon, lat in ll_inputs:
    x, y = tools.ll2gp(float(lat), float(lon), float(MODEL_LAT), float(MODEL_LON), float(MODEL_ROT), int(nx), int(ny), float(hh), float(dx_ts), float(dy_ts))
    print "%s %s %s %s %s" %(name, lon, lat, x, y)
    statcoordsfile.write('%5d %5d %5d %s\n'%(x,y,1,name))

llfile.close()
statcoordsfile.close()

#geoNet
```
Follow the examples listed in examples directory to download and process geoNet data.
Some of the products/plots produced during real-time ground motion simulation drills 
are dependent upon 
    Bradley_2010_Sa.py  
    calculateGMPE.py  
    readStationFile.py
which are copied over and placed in gmpe directory. 
The scripts that produce products for real-time simulation drills are placed 
in the examples directory. The example scripts make heavy use of utils.py and 
putils.py which contain many convenience functions. Note that some of the products
produced for real-time simulation drills are similar to those produced in 
groundMotionStationAnalysis code (now called post-processing).
```

```
Installation:
```
(0) copy setup.py so that geoNet and setup.py reside in the same parent directory then

(1) python2.7 setup.py build

(2) Separately compile rspectra.pyx  and manually place rpsectra.so in build/lib/geoNet
    Note the rspectra.so placed here has been compiled on hypocentre.

(3) python2.7 setup.py install --user 

  On hypocentre this places the geoNet package in ~/.local/lib64/python2.7/site-packages

```
Uninstall:
```
  rm -r ~/.local/lib64/python2.7/site-packages/geoNet

```
Usage example:
```
```python
from geoNet.utils import read_statsll

loc="/nesi/projects/nesi00213/RealTime/Obs/Mw5pt95_20150105_174841"
stats_fname="20150105_175001_eventStats_2016-12-12.ll"
  
event_stats = read_statsll(loc, stats_fname)
for stat_code in event_stats:
    lon, lat = event_stats[stat_code]
    print("{:^10s} {:^15.4f} {:^15.4f}".format(stats_code, lon, lat))

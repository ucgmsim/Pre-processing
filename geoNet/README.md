#geoNet
```
Follow the examples listed in examples directory to download and process geoNet data.
Some of the products/plots produced during real-time ground motion simulation drills 
The scripts that produce products for real-time simulation drills are placed 
in the examples directory. The example scripts make heavy use of utils.py and 
putils.py which contain many convenience functions. Note that some of the products
produced for real-time simulation drills are similar to those produced in 
groundMotionStationAnalysis code (now called post-processing).
```

```
Installation:

(1) python setup.py install --user 

  On hypocentre this places the geoNet package in ~/.local/lib64/python2.7/site-packages
```

```
Uninstall:

  rm -r ~/.local/lib64/python2.7/site-packages/geoNet
```

```
Usage example:
```
```python
from geoNet.geoNet.utils import read_statsll

loc="/nesi/projects/nesi00213/RealTime/Obs/Mw5pt95_20150105_174841"
stats_fname="20150105_175001_eventStats_2016-12-12.ll"
  
event_stats = read_statsll(loc, stats_fname)
for stat_code in event_stats:
    lon, lat = event_stats[stat_code]
    print("{:^10s} {:^15.4f} {:^15.4f}".format(stat_code, lon, lat))
```

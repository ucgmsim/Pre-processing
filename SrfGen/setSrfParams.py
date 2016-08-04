## Sets Variables for SRF/Stoch Generation

# TYPE:
# 1: point source to point source srf
# 2: point source to finite fault srf
# 3: finite fault to finite fault srf
TYPE = 2

# PREFIX for gsf/srf/stoch files
PREFIX = 'standard_'
# FILENAME for corners file
CORNERS = 'cnrs.txt'

###
### COMMON PARAMETERS (apply to all types)
###

# latitude (float)
LAT = -43.5029
# longitude (float)
LON = 172.8284
# depth (float)
DEPTH = 4.0
# magnitude (float)
MAG = 5.8
# strike (int)
STK = 54
# dip (int)
DIP = 75
# rake (int)
RAK = 137

###
### RELATING TO TYPE 1 (point to point)
###

# specify seismic moment directly (-1 to use magnitude)
MOM = -1

###
### RELATING TO TYPE 2,3 (finite fault output)
###
DT = 0.025
SEED = 1129571

###
### RELATING TO TYPE 2 (centroid moment tensor to finite fault)
###

# Mw Scaling Relation (string), one of:
# HanksBakun2002
# BerrymanEtAl2002
# VillamorEtAl2001
MWSR = 'BerrymanEtAl2002'

###
### RELATING TO TYPE 3 (ffd to finite fault)
###

FLEN = 16.0
DLEN = 0.10
FWID = 9.0
DWID = 0.10
DTOP = 0.63
# km east
SHYPO = -2.0
# km down
DHYPO = 6.0


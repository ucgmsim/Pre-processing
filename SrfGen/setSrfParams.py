## Sets Variables for SRF/Stoch Generation

# TYPE:
# 1: point source to point source srf
# 2: point source to finite fault srf
# 3: finite fault to finite fault srf
# 4: multi-segment finite fault srf
TYPE = 4

# PREFIX for gsf/srf/stoch files
PREFIX = 'standard_'
# FILENAME for corners file
CORNERS = 'cnrs.txt'

###
### COMMON PARAMETERS (apply to all types but multi)
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
### RELATING TO TYPES 2,3,4 (finite fault output)
###
DT = 0.025
SEED = 103245
# only used in type 4 to go over multiple seeds
SEED_INC = 1

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


###
### RELATING TO TYPE 4 (multi-segment finite fault)
###


# how many scenarios to run. different scenarios can have
# randomised parameters (randomisation amount set below)
N_SCENARIOS = 4
# how many scenarios to run before incrementing the seed value.
# this will randomise parameters but keep the seed the same.
N_SEED_INC = 2
# description of main segments
CASES = ['rv', 'ss']

# prefix for all files
M_NAME = 'bev01'
# master segments can be made up of multiple sub-segments
# set 
M_MAG = [6.50, 7.06]
# as above, if MOM invalid ( <= 0), calculated from MAG
M_MOM = [-1, -1]
# how many segments each master is made of
M_NSEG = [1, 6]

M_SEG_DELAY = ["0", "1"]
M_RVFAC_SEG = ["-1", "0.5,0.5"]
M_GWID = ["-1", "8.0,8.0"]
M_RUP_DELAY = ["0", "4.0"]

# 2 dimentional arrays for each set of segments
# demensions should match MASTER_NSEG
# if values are same, shortcut is [value] * repetitions
M_FLEN = [[10], [12, 20, 14, 7, 11, 8]]
M_DLEN = [[0.2], [0.2] * M_NSEG[1]]
M_FWID = [[18], [18] * M_NSEG[1]]
M_DWID = [[0.2], [0.2] * M_NSEG[1]]
M_DTOP = [[1.0], [0.0] * M_NSEG[1]]
M_STK = [[40], [121, 87, 87, 216, 40, 150]]
M_RAK = [[90], [-180, 180, 180, 78, 171, 6]]
M_DIP = [[75], [105, 85, 85, 50, 80, 54]]
M_ELON = [[172.135], [172.005, 172.195, 172.380, 171.9434, 172.3115, 172.1811]]
M_ELAT = [[-43.552], [-43.567, -43.590, -43.573, -43.5779, -43.5509, -43.5095]]
M_SHYPO = [[2], [6, -10, -27, 9.5, -25.5, 7]]
M_DHYPO = [[10], [6] * M_NSEG[1]]

### set variability, it will change values between scenarios
### length of arrays is same as CASES
# magnitude variability (absolute range)
V_MAG = [0.2, 0.5]
# fault width/length variability (proportion, eg. 0.1 = 10%)
V_FWID = [0.03, 0.15]
V_FLEN = [0.10, 0.00]

# TODO: magnitude as a function result if wanted

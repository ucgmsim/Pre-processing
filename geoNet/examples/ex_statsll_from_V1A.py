from pyqc import utils
import datetime
import os

loc="/nesi/projects/nesi00213/ObservedGroundMotions/ahsan/Mw4pt9_20110429_190804_11Jan2017/Vol1/data"
stats_fname="_".join(["20110429_190804", "eventStats", str(datetime.date.today())])
stats_fname+=".ll"
#utils.statsll_from_V1A(loc, save_stats=True)


loc_all_geoNet_stats="/nesi/projects/nesi00213/StationInfo"
fname_all_geoNet_stats="all_geoNet_stats+2016-12-20.ll"
loc_V1A=loc
utils.statsll_from_V1A_v2(loc_all_geoNet_stats, fname_all_geoNet_stats,
                         loc_V1A, save_stats=True,
                         fname=stats_fname, loc=os.getcwd())

import os
script_dir = os.path.dirname(os.path.abspath(__file__))

GSF_BIN='/nesi/projects/nesi00213/tools/fault_seg2gsf_dipdir'
# note version is missing, will be appended
FF_SRF_BIN='/nesi/projects/nesi00213/tools/genslip'
PS_SRF_BIN='/nesi/projects/nesi00213/tools/generic_slip2srf'
STOCH_BIN='/nesi/projects/nesi00213/tools/srf2stoch'
# default 1d velocity model
VELOCITY_MODEL='%s/lp_generic1d-gp01_v1.vmod' % (script_dir)


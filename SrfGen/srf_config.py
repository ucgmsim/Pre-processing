import os
script_dir = os.path.dirname(os.path.abspath(__file__))

# kupe
tools = '/nesi/transit/nesi00213/tools'

GSF_BIN='%s/fault_seg2gsf_dipdir' % (tools)
# note version is missing, will be appended
FF_SRF_BIN='%s/genslip' % (tools)
PS_SRF_BIN='%s/generic_slip2srf' % (tools)
STOCH_BIN='%s/srf2stoch' % (tools)
# default 1d velocity model
VELOCITY_MODEL='%s/lp_generic1d-gp01_v1.vmod' % (script_dir)

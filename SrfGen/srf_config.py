import os

script_dir = os.path.dirname(os.path.abspath(__file__))

# default 1d velocity model
VELOCITY_MODEL = "%s/lp_generic1d-gp01_v1.vmod" % (script_dir)

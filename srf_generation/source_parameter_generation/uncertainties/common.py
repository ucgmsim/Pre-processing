from collections import namedtuple
from typing import Any, Dict

GCMT_PARAM_NAMES = [
    "pid",  # str
    "lat",  # float
    "lon",  # float
    "depth",  # float
    "mag",  # float
    "strike",  # float
    "dip",  # float
    "rake",  # float
]
GCMT_Source = namedtuple("GCMT_Source", GCMT_PARAM_NAMES)

NHM_PARAM_NAMES = [
    "name",  # str
    "tect_type",  # str
    "fault_type",  # str
    "length_mean",  # float
    "length_sigma",  # float
    "dip_mean",  # float
    "dip_sigma",  # float
    "dip_dir",  # float
    "rake",  # float
    "dbot_mean",  # float
    "dbot_sigma",  # float
    "dtop_mean",  # float
    "dtop_min",  # float
    "dtop_max",  # float
    "slip_rate_mean",  # float
    "slip_rate_sigma",  # float
    "coupling_coefficient",  # float
    "coupling_coefficient_sigma",  # float
    "magnitude_mean",  # float
    "recurance_interval_mean",  # float
    "coordinates",  # List of (lon, lat) tuples
]
NHM_Source = namedtuple("NHM_Source", NHM_PARAM_NAMES)

GENERAL_PARAMS = ["type", "name"]

SRFGEN_TYPE_1_PARAMS = [
    "latitude",
    "longitude",
    "depth",
    "magnitude",
    "moment",
    "strike",
    "rake",
    "dip",
    "vs",
    "rho",
    "risetime",
    "dt",
    "stype",
    "inittime",
    "target_slip_cm",
    "target_area_km",
]

HF_RUN_PARAMS = [
    "sdrop",
    "t-sec",
    "rayset",
    "seed",
    "duration",
    "dt",
    "fmax",
    "kappa",
    "qfexp",
    "rvfac",
    "rvfac_shallow",
    "rvfac_deep",
    "czero",
    "calpha",
    "mom",
    "rupv",
    "velocity-model",
    "site-vm-dir",
    "vs-moho",
    "fa_sig1",
    "fa_sig2",
    "rv_sig1",
    "hf_path_dur",
]

LF_RUN_PARAMS = []

BB_RUN_PARAMS = ["flo", "fmin", "fmidbot", "lfvsref"]

RUN_TIME_PARAMS = HF_RUN_PARAMS + LF_RUN_PARAMS + BB_RUN_PARAMS


def verify_realisation_params(params: Dict[str, Any]):
    if params["type"] == 1:
        mismatch = [
            name
            for name in params.keys()
            if name not in GENERAL_PARAMS + SRFGEN_TYPE_1_PARAMS + RUN_TIME_PARAMS
        ]
    else:
        raise ValueError(
            f"'type' parameter given not valid. Given value {params['type']} is of type {type(params['type'])}."
        )
    if mismatch:
        raise ValueError(f"Unexpected parameters found: {mismatch}")

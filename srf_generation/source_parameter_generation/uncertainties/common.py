from collections import namedtuple
from typing import Any, Dict

GCMT_PARAM_NAMES = ["pid", "lat", "lon", "depth", "mag", "strike", "dip", "rake"]
GCMT_PARAMS = namedtuple("GCMT_Source", GCMT_PARAM_NAMES)


SRFGEN_TYPE_1_PARAMS = [
    "name",
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
    "rise_time",
    "rvfac",
    "rough",
    "slip_cov",
    "flen",
    "dlen",
    "fwid",
    "dwid",
    "dtop",
    "shypo",
    "dhypo",
]


def verify_params(params: Dict[str, Any]):
    mismatch = [name for name in params.keys() if name not in SRFGEN_TYPE_1_PARAMS]
    if mismatch:
        raise ValueError(f"Unexpected parameters found: {mismatch}")

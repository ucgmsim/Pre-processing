from typing import Any, Dict

SRFGEN_PARAMS = [
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
    mismatch = [name for name in params.keys() if name not in SRFGEN_PARAMS]
    if mismatch:
        raise ValueError(f"Unexpected parameters found: {mismatch}")

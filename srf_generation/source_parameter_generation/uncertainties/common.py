from collections import namedtuple
from typing import Any, Dict
from pandas import DataFrame
import numpy as np
from scipy.stats import randint

from srf_generation.pre_processing_common import calculate_corners, get_hypocentre
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    mw_2_lw_scaling_relation,
)

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

GENERAL_PARAMS = ["name", "type", "genslip_version", "srfgen_seed"]

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

SRFGEN_TYPE_2_PARAMS = [
    "latitude",
    "longitude",
    "depth",
    "magnitude",
    "moment",
    "strike",
    "rake",
    "dip",
    "rough",
    "slip_cov",
    "flen",
    "dlen",
    "fwid",
    "dwid",
    "dtop",
    "shypo",
    "dhypo",
    "rvfac",
    "mwsr",
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


def get_seed():
    """Returns a seed in the range of 0 to the largest 4 byte signed int possible in C"""
    return randint(0, 2 ** 31 - 1).rvs()


def filter_realisation_input_params(fault_type: int, params: Dict[str, Any]):
    if fault_type == 1:
        params = {
            key: value
            for key, value in params.items()
            if key in GENERAL_PARAMS + SRFGEN_TYPE_1_PARAMS + RUN_TIME_PARAMS
        }
    elif fault_type == 2:
        params = {
            key: value
            for key, value in params.items()
            if key in GENERAL_PARAMS + SRFGEN_TYPE_2_PARAMS + RUN_TIME_PARAMS
        }
    else:
        raise ValueError(
            f"'type' parameter given not valid. Given value {params['type']} is of type {type(params['type'])}."
        )
    return params


def verify_realisation_params(params: Dict[str, Any]):
    if params["type"] == 1:
        mismatch = [
            name
            for name in params.keys()
            if name not in GENERAL_PARAMS + SRFGEN_TYPE_1_PARAMS + RUN_TIME_PARAMS
        ]
    elif params["type"] == 2:
        mismatch = [
            name
            for name in params.keys()
            if name not in GENERAL_PARAMS + SRFGEN_TYPE_2_PARAMS + RUN_TIME_PARAMS
        ]
    else:
        raise ValueError(
            f"'type' parameter given not valid. Given value {params['type']} is of type {type(params['type'])}."
        )
    if mismatch:
        raise ValueError(f"Unexpected parameters found: {mismatch}")


def verify_1d_vel_mod(vel_mod_1d: DataFrame):
    assert vel_mod_1d.columns == ("depth", "vp", "vs", "rho", "qp", "qs")


def focal_mechanism_2_finite_fault(lat, lon, depth, mag, strike, rake, dip, mwsr):
    """
    Purpose: To create a finite fault geometry based on centroid moment tensor
    solution information and magnitude scaling relationships.
    The use of such finite fault geometry is for first order computation of
    source-to-site distances.

    revisions
    v2 - reoriented fault plane to be consistent with along strike direction
    v3 - added 'argout_emod3d' output and associated commands to output details
         which are input in the finite fault model generation for EMOD3D graves methodology

    assumptions:
    1) The centroid moment tensor represents the centroid of the fault plane,
    located half way along strike and downdip
    2) The fault area is computed from the median of a Moment scaling
    relationship and thus is assumed to be a 'typical' stress drop event for
    the Mw-relationship used.
    Should only be used for small events where the downdip width is smaller
    than the seismogenic depth (i.e. so the assumption of a square fault plane
    is reasonable); depth is reset to zero if above ground values determined
    3) srike in deg measured clockwise from N; rake in deg measured
    anti-clockwise from the strike direction (and in the range
    -180<lambda<180; dip in deg measured downward from horizontal

    Input variables:
    Lat    - the latitude of the focal mechanism (-ve below equator)
    Lon    - the lon of the focal mech (-180<lon<180)
    Depth  - the depth of the focal mech (in km)
    Mw     - the moment magnitude from the Mw tensor soln
    strike - the strike of the focal mech (only 1)
    rake   - rake (not used, but carried forward)
    dip    - the dip angle of the fault plane

    output variables:
    lat       - the latitude of the subfault
    lon       - the lon of the subfault
    depth     - depth of the subfault
    lonAsList - as a 1D array
    """

    # fixed values
    dlen = 0.1
    dwid = 0.1
    shypo = 0.00

    # get the fault geometry (square edge length)
    fault_length, fault_width = mw_2_lw_scaling_relation(mag, mwsr, rake)

    # number of subfaults
    nx = int(round(fault_length / dlen))
    ny = int(round(fault_width / dwid))

    # rounded subfault spacing
    dlen = fault_length / float(nx)
    dwid = fault_width / float(ny)

    # use cartesian coordinate system to define the along strike and downdip
    # locations taking the center of the fault plane as (x,y)=(0,0)
    x_pos: np.ndarray = np.arange(dlen / 2.0, fault_length, dlen) - fault_length / 2.0
    y_pos: np.ndarray = np.arange(dwid / 2.0, fault_width, dwid)[
        ::-1
    ] - fault_width / 2.0

    lats, lons = calculate_corners(dip, x_pos, y_pos, lat, lon, strike)

    depth_loc_relative = (-y_pos * np.sin(np.radians(dip))).repeat(ny).reshape((nx, ny))
    depths = np.maximum(depth + depth_loc_relative, 0)
    dhypo = fault_width / 2.0

    lon_tcl, lat_tcl = get_hypocentre(lat, lon, 0, -dhypo, strike, dip)

    depth_tcl = max([depth - dhypo * np.sin(np.radians(dip)), 0])

    return (
        lats,
        lons,
        depths,
        fault_length,
        dlen,
        fault_width,
        dwid,
        depth_tcl,
        lat_tcl,
        lon_tcl,
        shypo,
        dhypo,
    )

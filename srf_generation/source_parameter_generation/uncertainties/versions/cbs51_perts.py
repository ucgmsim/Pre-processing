"""Version 203 for Sarah Neills source parameter perturbations"""
from typing import Any, Dict

import numpy as np
import pandas as pd

from qcore import geo
from qcore.uncertainties import distributions
from qcore.uncertainties.mag_scaling import MagnitudeScalingRelations
from srf_generation.Fault import Type2
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
    get_seed,
    filter_realisation_input_params,
)

TYPE = 2

def uniform_dist(u_mean, u_half_range):
    """function for sampling uniform distribution"""
    u_sample = np.random.uniform(
        low=u_mean - u_half_range, high=u_mean + u_half_range, size=1
    )
    return u_sample

def generate_source_params(
    source_data: GCMT_Source,
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame,
    vs30_data: pd.DataFrame = None,
    **kwargs,
) -> Dict[str, Any]:
    additional_source_parameters = filter_realisation_input_params(
        TYPE, additional_source_parameters
    )

    realisation = kwargs

    params = generate_from_gcmt(source_data)
    params.update(additional_source_parameters)
    params.update({"dt": 0.005, "seed": get_seed(), "genslip_version": "5.4.2"})
    realisation["params"] = params

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    verify_realisation_params(params)
    return realisation


def generate_from_gcmt(source_data: GCMT_Source):

    mag = distributions.truncated_normal(source_data.mag, 0.075, 2)  # magnitude
    lat_temp, lon_temp = geo.ll_shift(
        source_data.lat, source_data.lon, distributions.truncated_normal(0.0, 1.0, 2), 0
    )
    lat, lon = geo.ll_shift(
        lat_temp, lon_temp, distributions.truncated_normal(0.0, 1.0, 2), 90
    )
    depth = distributions.truncated_normal(source_data.depth, 2.0, 2)
    strike = distributions.truncated_normal(source_data.strike, 10, 2)
    dip = distributions.truncated_normal(source_data.dip, 10, 2)
    rake = distributions.truncated_normal(source_data.rake, 15, 4)

    # correct strike dip and rake:
    if dip > 90:
        strike = (strike + 180) % 360  # flip the strike
        dip = 180 - dip
        rake = rake * -1

    if dip < 0:
        strike = (strike + 180) % 360
        dip = 90 - dip % 90
        rake = rake * -1

    if rake > 180:
        rake = rake - 360

    if rake < -180:
        rake = rake + 360

    if strike > 360:
        strike = strike % 360

    if strike < 0:
        strike = strike % 360

    ### the following parameters feed into the sim_params.yaml

    sdrop = distributions.truncated_log_normal(50, 0.3, 2)
    rvfac = uniform_dist(0.8, 0.075)  # rupture velocity factor
    kappa = distributions.truncated_log_normal(0.045, 0.3, 2)

    fault: Type2 = Type2(
        source_data.pid,
        lat,
        lon,
        mag,
        strike,
        rake,
        dip,
        depth,
    )
    fault.magnitude_scaling_relation = MagnitudeScalingRelations.LEONARD2014

    params = {
        "kappa": kappa,
        "sdrop": sdrop,
        "rvfac": rvfac,
    }
    params.update(fault.to_dict())
    return params

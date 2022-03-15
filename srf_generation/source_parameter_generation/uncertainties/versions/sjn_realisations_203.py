"""Version 203 for Sarah Neills source parameter perturbations"""
from typing import Any, Dict

import numpy as np

from qcore import geo
from qcore.uncertainties import distributions
from qcore.uncertainties.mag_scaling import mag2mom
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
)


def generate_source_params(source_data: GCMT_Source) -> Dict[str, Any]:
    """source_data should have the following parameters available via . notation:
    - source_data.pid: name of the event
    - source_data.lat: latitude
    - source_data.lon: longitude
    - source_data.depth: depth
    - source_data.mag: magnitude
    - source_data.strike
    - source_data.dip
    - source_data.rake
    """

    params = generate_from_gcmt(source_data)

    verify_realisation_params(params)
    return params


# distribution functions:


def uniform_dist(u_mean, u_half_range):
    """function for sampling uniform distribution"""
    u_sample = np.random.uniform(
        low=u_mean - u_half_range, high=u_mean + u_half_range, size=1
    )
    return u_sample


def generate_from_gcmt(source_data: GCMT_Source):

    # area = mw_2_a_scaling_relation(
    #    source_data.mag,
    #    MagnitudeScalingRelations.LEONARD2014.value,
    #    source_data.strike,
    # )

    ### the following parameters feed into srfgen

    mag = distributions.truncated_normal(source_data.mag, 0.075, 2)  # magnitude
    # rvfrac = uniform_dist(0.8, 0.075)                   #rupture velocity factor
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

    # calculate normalised perturbations:
    z_mw_pert = (mag - source_data.mag) / 0.075
    z_lat_pert = (
        geo.ll_dist(source_data.lon, source_data.lat, source_data.lon, lat) / 1.0
    )
    z_long_pert = (
        geo.ll_dist(source_data.lon, source_data.lat, lon, source_data.lat) / 1.0
    )
    z_depth_pert = (depth - source_data.depth) / 2.0
    z_strike_pert = (strike - source_data.strike) / 10.0
    z_dip_pert = (dip - source_data.dip) / 10.0
    z_rake_pert = (rake - source_data.rake) / 15.0

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

    # calculate normalised perturbations:
    z_sdrop_pert = (sdrop - 50) / 1.35
    z_rvfac_pert = (rvfac - 0.8) / 0.075
    z_kappa_pert = (kappa - 0.045) / 1.35

    params = {
        "type": 1,
        "name": source_data.pid,
        "latitude": lat,
        "longitude": lon,
        "depth": depth,
        "magnitude": mag,
        "moment": float(mag2mom(mag)),
        "strike": strike,
        "dip": dip,
        "rake": rake,
        "kappa": kappa,
        "sdrop": sdrop,
        "rvfac": rvfac,
    }
    return params

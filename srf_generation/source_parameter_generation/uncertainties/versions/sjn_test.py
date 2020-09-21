"""The template for future perturbation versions.
Update this docstring with information about the version"""
import numpy as np
import pandas as pd

from srf_generation.source_parameter_generation.uncertainties import distributions
from qcore import geo
from typing import Any, Dict, Union  # new line

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    mag2mom,
)


def verify_1d_vel_mod(vel_mod_1d: pd.DataFrame):  # new line
    for line in vel_mod_1d.itertuples():  # new line
        pass
        # assert line._fields == ('depth', 'vp', 'vs', 'rho', 'qp', 'qs')       #new line


def generate_source_params(
    source_data: Union[GCMT_Source],
    additional_source_parameters: Dict[str, Any],  # new line
    vel_mod_1d: pd.DataFrame = None,  # new line
    hf_vel_mod_1d: pd.DataFrame = None,  # new line
    vs30_data: pd.DataFrame = None,
    **kwargs,  # new line
) -> Dict[str, Any]:
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

    realisation = kwargs  # new line
    ### Start of custom code area

    if isinstance(source_data, GCMT_Source):
        realisation["params"] = generate_from_gcmt(
            source_data, additional_source_parameters
        )  # modified
    else:
        raise ValueError("Not a valid type")

    if vel_mod_1d is not None:  # new line
        vs_profile = distributions.truncated_log_normal(vel_mod_1d.vs, 0.05, 4)
        vel_mod_1d.vs = vs_profile

        realisation["vel_mod_1d"] = vel_mod_1d

    # HF 1d profile
    if hf_vel_mod_1d is not None or True:
        # HF_onedprofile_pert = vel_mod_1d
        hf_vel_mod_1d = pd.read_csv(
            "/home/jpa198/tmp/Cant1D_v3-midQ_OneRay.1d",
            delim_whitespace=True,
            names=["depth", "vp", "vs", "rho", "qp", "qs"],
            skiprows=1,
        )
        vals = distributions.truncated_log_normal(hf_vel_mod_1d.qs, 0.3, 2.5)
        hf_vel_mod_1d.qs = vals

        realisation["hf_vel_mod_1d"] = hf_vel_mod_1d

    if vs30_data is not None:
        # print(f"Got vs30: {vs30_data}")
        realisation["vs30"] = vs30_data.copy(deep=True)
        realisation["vs30"]["vs30"] = distributions.truncated_log_normal(
            vs30_data["median"].values, vs30_data["sigma"].values, 2
        )
    # else:
    # print("Didn't get vs30")

    ### End of custom code area

    verify_realisation_params(realisation["params"])  # modified

    if vel_mod_1d is not None:  # new
        verify_1d_vel_mod(realisation["vel_mod_1d"])  # new

    return realisation  # modified


# distribution functions:


def uniform_dist(u_mean, u_half_range):
    """function for sampling uniform distribution"""
    u_sample = np.random.uniform(
        low=u_mean - u_half_range, high=u_mean + u_half_range, size=1
    )
    return u_sample


def generate_from_gcmt(
    sources_line: GCMT_Source, additional_source_parameters: Dict[str, Any]
):  # modified

    # area = mw_2_a_scaling_relation(
    #    sources_line.mag,
    #    MagnitudeScalingRelations.LEONARD2014.value,
    #    sources_line.strike,
    # )

    ### the following parameters feed into srfgen

    mag = distributions.truncated_normal(sources_line.mag, 0.075, 2)  # magnitude
    # rvfrac = uniform_dist(0.8, 0.075)                   #rupture velocity factor
    lat_temp, lon_temp = geo.ll_shift(
        sources_line.lat,
        sources_line.lon,
        distributions.truncated_normal(0.0, 1.0, 2),
        0,
    )
    lat, lon = geo.ll_shift(
        lat_temp, lon_temp, distributions.truncated_normal(0.0, 1.0, 2), 90
    )
    depth = distributions.truncated_normal(sources_line.depth, 2.0, 2)
    strike = distributions.truncated_normal(sources_line.strike, 10, 2)
    dip = distributions.truncated_normal(sources_line.dip, 10, 2)
    rake = distributions.truncated_normal(sources_line.rake, 15, 4)
    qsfrac = distributions.truncated_log_normal(50, 0.3, 2.5)

    # calculate normalised perturbations:
    z_mw_pert = (mag - sources_line.mag) / 0.075
    z_lat_pert = (
        geo.ll_dist(sources_line.lon, sources_line.lat, sources_line.lon, lat) / 1.0
    )
    z_long_pert = (
        geo.ll_dist(sources_line.lon, sources_line.lat, lon, sources_line.lat) / 1.0
    )
    z_depth_pert = (depth - sources_line.depth) / 2.0
    z_strike_pert = (strike - sources_line.strike) / 10.0
    z_dip_pert = (dip - sources_line.dip) / 10.0
    z_rake_pert = (rake - sources_line.rake) / 15.0

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
        "name": sources_line.pid,
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
        "qsfrac": qsfrac,
    }

    params.update(additional_source_parameters)  # modified
    return params

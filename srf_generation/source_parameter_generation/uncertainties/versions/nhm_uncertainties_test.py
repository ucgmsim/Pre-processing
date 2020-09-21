#!/usr/bin/env python
"""The template for future perturbation versions.
Update this docstring with information about the version"""
import pandas as pd
from qcore.nhm import NHMFault

from srf_generation.Fault import Type4, fault_factory
from srf_generation.source_parameter_generation.uncertainties import distributions
from typing import Any, Dict

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    get_seed,
)


def generate_source_params(
    source_data: NHMFault,
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame = None,
    vs30_data: pd.DataFrame = None,
    **kwargs,
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

    realisation = kwargs

    realisation["params"] = generate_params_from_nhm(
        source_data, additional_source_parameters
    )

    if vel_mod_1d is not None:
        onedprofile_pert = vel_mod_1d.copy()
        vs_profile = [
            distributions.truncated_log_normal(vs, 0.05, 4)
            for vs in onedprofile_pert.vs
        ]
        onedprofile_pert.vs = vs_profile

        realisation["vel_mod_1d"] = onedprofile_pert

    # HF 1d profile
    HF_onedprofile_ori = pd.read_csv(
        "/nesi/project/nesi00213/VelocityModel/Mod-1D/Cant1D_v3-midQ_OneRay.1d",
        delim_whitespace=True,
        names=["depth", "vp", "vs", "rho", "qp", "qs"],
        skiprows=1,
    )
    HF_onedprofile_pert = HF_onedprofile_ori.copy()
    HF_qs_profile = [
        distributions.truncated_log_normal(qs, 0.3, 2.5)
        for qs in HF_onedprofile_pert.qs
    ]
    HF_onedprofile_pert.qs = HF_qs_profile

    realisation["hf_vel_mod_1d"] = HF_onedprofile_pert

    if vs30_data is not None:
        realisation["vs30"] = vs30_data.copy(deep=True)
        realisation["vs30"]["vs30"] = distributions.truncated_log_normal(
            vs30_data["median"].values, vs30_data["sigma"].values, 2
        )

    return realisation


def generate_params_from_nhm(
    source_data: NHMFault, additional_source_parameters: Dict[str, Any]
):

    fault: Type4 = fault_factory(4)(source_data)

    fault.shypo = fault.length * distributions.rand_shyp()
    fault.dhypo = fault.width * distributions.truncated_weibull(1)

    # the following parameters feed into srfgen
    # Start changing geometry
    fault.magnitude = distributions.truncated_normal(
        fault.magnitude, fault.magnitude_sigma, 2
    )

    fault.dip = distributions.truncated_normal(fault.dip, source_data.dip_sigma, 2)
    fault.dip_dir = distributions.truncated_normal(fault.dip_dir, 10, 2)
    fault.rake = distributions.truncated_normal(fault.rake, 15, 4)

    # the following parameters feed into the sim_params.yaml
    sdrop = distributions.truncated_log_normal(50, 0.3, 2)
    rvfac = distributions.uniform(0.8, 0.075)  # rupture velocity factor
    kappa = distributions.truncated_log_normal(0.045, 0.3, 2)

    params = fault.to_dict()
    params.update(
        {
            "dt": 0.005,
            "seed": get_seed(),
            "genslip_version": "5.4.2",
            "sdrop": sdrop,
            "rvfac": rvfac,
            "kappa": kappa,
        }
    )

    params.update(additional_source_parameters)

    verify_realisation_params(params)

    return params

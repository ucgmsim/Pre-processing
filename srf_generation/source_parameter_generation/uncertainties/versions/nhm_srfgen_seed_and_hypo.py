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

    return realisation


def generate_params_from_nhm(
    source_data: NHMFault, additional_source_parameters: Dict[str, Any]
):

    fault: Type4 = fault_factory(4)(source_data)

    # Always in the middle of the fault
    fault.shypo = fault.length * distributions.rand_shyp()
    fault.dhypo = fault.width * distributions.truncated_weibull(1)

    params = fault.to_dict()
    params.update(
        {
            "dt": 0.005,
            "seed": get_seed(),
            "genslip_version": "5.4.2",
        }
    )

    params.update(additional_source_parameters)

    verify_realisation_params(params)

    return params

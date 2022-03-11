"""A basic perturbator as an example and starting point"""
from typing import Any, Dict

import numpy as np
import pandas as pd

from qcore.nhm import NHMFault
from qcore.uncertainties.mag_scaling import lw_to_mw_sigma_scaling_relation
from srf_generation.Fault import fault_factory, Type4
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    get_seed,
)
from srf_generation.source_parameter_generation.uncertainties.distributions import (
    rand_shyp,
    truncated_normal,
    truncated_weibull,
)

TYPE = 4


def generate_source_params(
    source_data: NHMFault,
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame,
    vs30_data: pd.DataFrame = None,
    **kwargs
) -> Dict[str, Any]:
    """source_data should have the following parameters available via . notation:
    - source_data.pid: name of the event
    - source_data.lat: latitude
    - source_data.lon: longitude
    - source_data.depth
    - source_data.mag: magnitude
    - source_data.strike
    - source_data.dip
    - source_data.rake
    """

    fault: Type4 = fault_factory(TYPE)(source_data)

    fault.shypo = (fault.length / 2) * rand_shyp()
    fault.dhypo = fault.width * truncated_weibull(1)

    fault.rake = truncated_normal(fault.rake, 15, 4)
    mag, sigma = lw_to_mw_sigma_scaling_relation(
        fault.length, fault.width, fault.mwsr, fault.rake
    )
    perturbated_magnitude = truncated_normal(mag, sigma, 1)

    params = fault.to_dict()
    params.update({"dt": 0.005, "seed": get_seed(), "genslip_version": "5.4.2"})

    params["sdrop"] = 50 * np.sqrt(10 ** (perturbated_magnitude - mag))
    params["magnitude"] = perturbated_magnitude

    realisation = kwargs

    realisation["params"] = params
    realisation["params"].update(additional_source_parameters)

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    verify_realisation_params(realisation["params"])
    return realisation

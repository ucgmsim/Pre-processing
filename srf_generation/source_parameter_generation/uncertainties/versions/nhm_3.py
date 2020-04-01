"""A basic perturbator as an example and starting point"""
import pandas as pd
from typing import Any, Dict

from qcore.nhm import NHMFault

from srf_generation.Fault import fault_factory, Type4
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
)

TYPE = 3


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

    params = fault.to_dict()

    realisation = kwargs

    realisation["params"] = params
    realisation["params"].update(additional_source_parameters)

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    verify_realisation_params(realisation["params"])
    return realisation

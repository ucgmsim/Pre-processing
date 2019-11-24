"""A basic perturbator as an example and starting point"""
import pandas as pd
from typing import Any, Dict, Union

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
    NHM_Source,
)


def generate_source_params(
    source_data: Union[GCMT_Source, NHM_Source],
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame,
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

    realisation = kwargs

    realisation["params"] = {
        "type": 1,
        "name": source_data.pid,
        "latitude": source_data.lat,
        "longitude": source_data.lon,
        "depth": source_data.depth,
        "magnitude": source_data.mag,
        "strike": source_data.strike,
        "dip": source_data.dip,
        "rake": source_data.rake,
    }
    realisation["params"].update(additional_source_parameters)

    verify_realisation_params(realisation["params"])
    return realisation

"""A basic perturbator as an example and starting point"""
from typing import Any, Dict, Union

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_params,
    GCMT_Source,
    NHM_Source,
)


def generate_source_params(
    sources_line: Union[GCMT_Source, NHM_Source],
    additional_source_parameters: Dict[str, Any],
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

    params = {
        "type": 1,
        "name": sources_line.pid,
        "latitude": sources_line.lat,
        "longitude": sources_line.lon,
        "depth": sources_line.depth,
        "magnitude": sources_line.mag,
        "strike": sources_line.strike,
        "dip": sources_line.dip,
        "rake": sources_line.rake,
    }
    params.update(additional_source_parameters)

    verify_params(params)
    return params

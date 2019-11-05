"""A basic perturbator as an example and starting point"""
import numpy as np
from typing import Any, Dict, Union

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_params,
    GCMT_Source,
    NHM_Source,
)
from srf_generation.source_parameter_generation.uncertainties.distributions import (
    truncated_normal,
    truncated_log_normal,
    uniform,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    mw_2_a_scaling_relation,
    MagnitudeScalingRelations,
    mag2mom,
)


def generate_source_params(
    sources_line: Union[GCMT_Source, NHM_Source],
    additional_source_parameters: Dict[str, Any],
    **kwargs,
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

    magnitude = truncated_normal(mean=sources_line.mag, std_dev=0.05, std_dev_limit=2)

    params = {
        "type": 1,
        "name": sources_line.pid,
        "latitude": sources_line.lat,
        "longitude": sources_line.lon,
        "depth": sources_line.depth,
        "magnitude": magnitude,
        "moment": float(mag2mom(magnitude)),
        "strike": sources_line.strike,
        "dip": sources_line.dip,
        "rake": sources_line.rake,
        "sdrop": float(truncated_log_normal(mean=50, std_dev=0.3, std_dev_limit=2)),
        "risetime": float(uniform(mean=0.8, half_range=0.075)),
    }

    verify_params(params)
    return params

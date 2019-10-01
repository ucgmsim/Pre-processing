"""A basic perturbator as an example and starting point"""
import numpy as np
from typing import Any, Dict

from SrfGen.source_parameter_generation.uncertainties.common import (
    verify_params,
    GCMT_PARAMS,
)
from SrfGen.source_parameter_generation.uncertainties.mag_scaling import (
    mw_2_a_scaling_relation,
    MagnitudeScalingRelations,
    mag2mom,
)


def generate_source_params(sources_line: GCMT_PARAMS) -> Dict[str, Any]:
    """sources_line should have the following parameters available via . notation:
      - sources_line.pid: name of the event
      - sources_line.lat: latitude
      - sources_line.lon: longitude
      - sources_line.depth
      - sources_line.mag: magnitude
      - sources_line.strike
      - sources_line.dip
      - sources_line.rake
    """

    area = mw_2_a_scaling_relation(
        sources_line.mag,
        MagnitudeScalingRelations.LEONARD2014.value,
        sources_line.strike,
    )

    params = {
        "type": 1,
        "name": sources_line.pid,
        "latitude": sources_line.lat,
        "longitude": sources_line.lon,
        "depth": sources_line.depth,
        "magnitude": sources_line.mag,
        "moment": mag2mom(sources_line.mag),
        "strike": sources_line.strike,
        "dip": sources_line.dip,
        "rake": sources_line.rake,
        "width": np.sqrt(area),
        "length": np.sqrt(area),
    }

    verify_params(params)
    return params

"""A basic perturbator as an example and starting point"""
from typing import Any, Dict, Union

from pandas import DataFrame
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
    NHM_Source,
)
from srf_generation.source_parameter_generation.uncertainties.distributions import (
    truncated_normal,
    truncated_log_normal,
    uniform,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import mag2mom


def generate_source_params(
    sources_line: Union[GCMT_Source, NHM_Source],
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: DataFrame = None,
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

    realisation = kwargs

    magnitude = truncated_normal(mean=sources_line.mag, std_dev=0.05, std_dev_limit=2)

    realisation["params"] = {
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
    realisation["params"].update(additional_source_parameters)
    if vel_mod_1d is not None:
        realisation["vel_mod_1d"] = vel_mod_1d

    verify_realisation_params(realisation["params"])
    return realisation

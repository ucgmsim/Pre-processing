"""A basic perturbator as an example and starting point"""
from typing import Any, Dict

from pandas import DataFrame
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
)
from srf_generation.source_parameter_generation.uncertainties.distributions import (
    truncated_normal,
    truncated_log_normal,
    uniform,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import mag2mom


def generate_source_params(
    source_data: GCMT_Source,
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: DataFrame,
    vs30_data: DataFrame = None,
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

    magnitude = truncated_normal(mean=source_data.mag, std_dev=0.05, std_dev_limit=2)

    realisation["params"] = {
        "type": 1,
        "name": source_data.pid,
        "latitude": source_data.lat,
        "longitude": source_data.lon,
        "depth": source_data.depth,
        "magnitude": magnitude,
        "moment": float(mag2mom(magnitude)),
        "strike": source_data.strike,
        "dip": source_data.dip,
        "rake": source_data.rake,
        "sdrop": float(truncated_log_normal(mean=50, std_dev=0.3, std_dev_limit=2)),
        "risetime": float(uniform(mean=0.8, half_range=0.075)),
    }
    realisation["params"].update(additional_source_parameters)

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    verify_realisation_params(realisation["params"])
    return realisation

"""A basic perturbator as an example and starting point"""
import pandas as pd
from typing import Any, Dict

from srf_generation.Fault import fault_factory, Type1
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
)

TYPE = 1


def generate_source_params(
    source_data: GCMT_Source,
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

    fault: Type1 = fault_factory(TYPE)(
        source_data.pid,
        source_data.lat,
        source_data.lon,
        source_data.mag,
        source_data.strike,
        source_data.rake,
        source_data.dip,
        source_data.depth,
    )

    params = fault.to_dict()

    realisation = kwargs

    realisation["params"] = params
    realisation["params"].update(additional_source_parameters)

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]


    bad_params = verify_realisation_params(realisation["params"], throw_exception=False)
    for bad_param in bad_params:
        print("WARNING: Unsupported parameters will be ignored {bad_param}")
        realisation["params"].pop(bad_param)

    return realisation

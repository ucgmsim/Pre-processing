"""The template for future perturbation versions.
Update this docstring with information about the version"""
import pandas as pd
from typing import Any, Dict

from srf_generation.Fault import fault_factory, Type2
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
    get_seed,
    filter_realisation_input_params,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    MagnitudeScalingRelations,
)

TYPE = 2


def generate_source_params(
    source_data: GCMT_Source,
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame,
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
    additional_source_parameters = filter_realisation_input_params(
        TYPE, additional_source_parameters
    )

    realisation = kwargs

    fault: Type2 = fault_factory(TYPE)(
        source_data.pid,
        source_data.lat,
        source_data.lon,
        source_data.mag,
        source_data.strike,
        source_data.rake,
        source_data.dip,
        source_data.depth,
    )
    if "tect_type" in additional_source_parameters.keys() and additional_source_parameters["tect_type"] == "SUBDUCTION_INTERFACE":
        fault.magnitude_scaling_relation = MagnitudeScalingRelations.SKARLATOUDIS2016
    else:
        fault.magnitude_scaling_relation = MagnitudeScalingRelations.LEONARD2014

    params = fault.to_dict()

    params.update({"dt": 0.005, "seed": get_seed(), "genslip_version": "3.3"})

    params.update(additional_source_parameters)
    realisation["params"] = params

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    # End of custom code area
    verify_realisation_params(realisation["params"])
    return realisation

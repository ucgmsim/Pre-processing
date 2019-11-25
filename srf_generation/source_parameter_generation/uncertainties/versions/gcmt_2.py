"""The template for future perturbation versions.
Update this docstring with information about the version"""
import pandas as pd
from typing import Any, Dict, Union

from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
    NHM_Source,
    focal_mechanism_2_finite_fault,
    get_seed,
    filter_realisation_input_params,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    MagnitudeScalingRelations,
)

TYPE = 2


def generate_source_params(
    source_data: Union[GCMT_Source, NHM_Source],
    additional_source_parameters: Dict[str, Any],
    vel_mod_1d: pd.DataFrame = None,
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
    additional_source_parameters = filter_realisation_input_params(TYPE, additional_source_parameters)

    realisation = kwargs

    # Start of custom code area
    mwsr = MagnitudeScalingRelations.LEONARD2014
    if "mwsr" in additional_source_parameters.keys():
        mwsr = MagnitudeScalingRelations(additional_source_parameters.pop("mwsr"))

    flen, dlen, fwid, dwid, dtop, lat0, lon0, shypo, dhypo = focal_mechanism_2_finite_fault(
        source_data.lat,
        source_data.lon,
        source_data.depth,
        source_data.mag,
        source_data.strike,
        source_data.rake,
        source_data.dip,
        mwsr,
    )[
        3:
    ]

    params = {
        "type": TYPE,
        "name": source_data.pid,
        "latitude": source_data.lat,
        "longitude": source_data.lon,
        "depth": source_data.depth,
        "magnitude": source_data.mag,
        "strike": source_data.strike,
        "dip": source_data.dip,
        "rake": source_data.rake,
        "flen": flen,
        "dlen": dlen,
        "fwid": fwid,
        "dwid": dwid,
        "dtop": dtop,
        "shypo": shypo,
        "dhypo": dhypo,
        "dt": 0.005,
        "seed": get_seed(),
        "genslip_version": "3.3",
        "mwsr": mwsr.name,
    }

    params.update(additional_source_parameters)
    realisation["params"] = params
    # End of custom code area
    verify_realisation_params(realisation["params"])
    return realisation

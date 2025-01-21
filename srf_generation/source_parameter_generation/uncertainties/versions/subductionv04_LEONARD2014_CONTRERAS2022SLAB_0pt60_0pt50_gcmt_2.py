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

def _compute_dhypo(tectclass,mag,hypo_offset=float('0pt60'.replace('pt','.'))):

    if tectclass == "SUBDUCTION_INTERFACE":
        if mag <= 5:
            dhypo = 0.5
        elif mag <= 8:
            dhypo = 0.5 + (mag-5)*(hypo_offset-0.5)/3
        else:
            dhypo = hypo_offset
    else:
        dhypo = 0.5

    return dhypo

def _get_scaling_relation(tectclass,relations,interface='LEONARD2014',slab='CONTRERAS2022SLAB'):

    relation_dict = {
    
        # INTERFACE MOELS
        "LEONARD2014":relations.LEONARD2014,
        "SKARLATOUDIS2016":relations.SKARLATOUDIS2016,
        "CONTRERAS2022INTERFACE":relations.CONTRERAS2022INTERFACE,
        "TEST2022INTERFACE":relations.TEST2022INTERFACE,
        "THINGBAIJAM2017":relations.THINGBAIJAM2017,
        "BLASER2010":relations.BLASER2010,
        "MUROTANI2013":relations.MUROTANI2013,
        "ALLEN2017INTERFACELINEAR":relations.ALLEN2017INTERFACELINEAR,
        "ALLEN2017INTERFACEBILINEAR":relations.ALLEN2017INTERFACEBILINEAR,
        "STRASSER2010INTERFACE":relations.STRASSER2010INTERFACE,

        # SLAB MODELS
        "LEONARD2014":relations.LEONARD2014,
        "ALLEN2017SLAB":relations.ALLEN2017SLAB,
        "STRASSER2010SLAB":relations.STRASSER2010SLAB,
        "CONTRERAS2022SLAB":relations.CONTRERAS2022SLAB,
        "TEST2022SLAB":relations.TEST2022SLAB,
        }


    if tectclass == "SUBDUCTION_INTERFACE":
        out = relation_dict[interface]
    elif tectclass == "SUBDUCTION_SLAB":
        out = relation_dict[slab]
    else:
        out = relations.LEONARD2014

    return out

# _get_scaling_relation(additional_source_parameters["tect_type"],MagnitudeScalingRelations)

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
    if ("tect_type" in additional_source_parameters.keys() and additional_source_parameters["tect_type"] in ["SUBDUCTION_INTERFACE","SUBDUCTION_SLAB"]):
        fault.magnitude_scaling_relation = _get_scaling_relation(additional_source_parameters["tect_type"],MagnitudeScalingRelations)
        # fault.magnitude_scaling_relation = MagnitudeScalingRelations.SKARLATOUDIS2016
        #fault.magnitude_scaling_relation = MagnitudeScalingRelations.LEONARD2014 
    else:
        fault.magnitude_scaling_relation = MagnitudeScalingRelations.LEONARD2014

    params = fault.to_dict()

    params.update({"dt": 0.005, "seed": get_seed(), "genslip_version": "5.4.2"})

    if ("tect_type" in additional_source_parameters.keys() and additional_source_parameters["tect_type"] in ["SUBDUCTION_INTERFACE","SUBDUCTION_SLAB"]):
        params.update({"dhypo":_compute_dhypo(additional_source_parameters["tect_type"],params['magnitude'])}) 
    else:
        params.update({"dhypo":0.5})


    params.update(additional_source_parameters)
    realisation["params"] = params
    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    # End of custom code area
    verify_realisation_params(realisation["params"])
    return realisation

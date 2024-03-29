"""A basic perturbator as an example and starting point"""
import pandas as pd
from typing import Any, Dict

from qcore.nhm import NHMFault

from srf_generation import Asperity
from srf_generation.Fault import fault_factory, Type4
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    get_seed,
)
from srf_generation.source_parameter_generation.uncertainties.distributions import (
    rand_shyp,
    truncated_weibull,
)

TYPE = 4


def generate_source_params(
    source_data: NHMFault,
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

    fault: Type4 = fault_factory(TYPE)(source_data)

    fault.shypo = fault.length / 2 * rand_shyp()
    fault.dhypo = fault.width * truncated_weibull(1)

    params = fault.to_dict()
    params.update({"dt": 0.005, "seed": get_seed(), "genslip_version": "5.4.2"})

    left = -fault.length / 4
    right = fault.length / 4
    top = fault.width / 4
    bottom = 3 * fault.width / 4
    bottom_asperity = Asperity.Asperity(left, right, top, bottom, v=3)
    top = 3 * fault.width / 8
    bottom = 5 * fault.width / 8
    top_asperity = Asperity.Asperity(left / 2, right / 2, top, bottom, v=9)

    realisation = kwargs

    realisation["asperities"] = {
        "background": 1,
        "asperities": [bottom_asperity, top_asperity],
    }

    realisation["params"] = params
    realisation["params"].update(additional_source_parameters)

    if vs30_data is not None:
        realisation["vs30"] = vs30_data
        realisation["vs30"]["vs30"] = vs30_data["median"]

    verify_realisation_params(realisation["params"])
    return realisation

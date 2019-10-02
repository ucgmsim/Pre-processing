import pytest
import numpy as np

from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    mw_2_a_scaling_relation,
    a_2_mw_scaling_relation,
)

mwsrs = ["HanksBakun2002", "BerrymanEtAl2002", "VillamorEtAl2001", "Leonard2014"]


@pytest.mark.parametrize(
    ["expected_value", "rake"], [(0.1, 60), (2.0, 60), (3.0, 60), (4.0, 60), (5.0, 60)]
)
def test_mwsr_functions(expected_value, rake):
    for mwsr in mwsrs:
        assert np.isclose(
            a_2_mw_scaling_relation(
                mw_2_a_scaling_relation(expected_value, mwsr, rake), mwsr, rake
            ),
            expected_value,
        )

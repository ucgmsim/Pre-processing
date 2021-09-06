import pytest
import numpy as np

from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    mw_to_lw_scaling_relation,
    lw_to_mw_scaling_relation,
    MagnitudeScalingRelations,
)

mwsrs = MagnitudeScalingRelations


@pytest.mark.parametrize(
    ["expected_value", "rake"], [(0.1, 60), (2.0, 60), (3.0, 60), (4.0, 60), (5.0, 60)]
)
def test_mwsr_functions(expected_value, rake):
    for mwsr in mwsrs:
        assert np.isclose(
            lw_to_mw_scaling_relation(
                *mw_to_lw_scaling_relation(expected_value, mwsr, rake), mwsr, rake
            ),
            expected_value,
        )

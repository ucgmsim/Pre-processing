from enum import Enum
from typing import Union

import numpy as np


class MagnitudeScalingRelations(Enum):
    HANKSBAKUM2002 = "HanksBakun2002"
    BERRYMANETAL2002 = "BerrymanEtAl2002"
    VILLAMORETAL2001 = "VillamorEtAl2001"
    LEONARD2014 = "Leonard2014"


def mw_2_a_scaling_relation(mw: float, mw_scaling_rel: str, rake: Union[None, float] = None, leonard_ds: float = 4.00, leonard_ss: float = 3.99):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUM2002.value:
        try:
            # only small magnitude case
            assert mw <= 6.71
        except AssertionError as e:
            e.args += ("Cannot use HanksAndBakun2002 equation for mw > 6.71",)
            raise
        A = 10 ** (mw - 3.98)

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002.value:
        A = 10 ** (mw - 4.18)
        # i.e. set L=W=sqrt(A) in Eqn 1 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001.value:
        A = 10 ** (0.75 * (mw - 3.39))
        # i.e. Eqn 2 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014.value:
        if round(rake % 360 / 90.0) % 2:
            A = 10 ** (mw - leonard_ds)
        else:
            A = 10 ** (mw - leonard_ss)

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # Area
    return A


def a_2_mw_scaling_relation(a, mw_scaling_rel, rake=None, leonard_ds=4.00, leonard_ss=3.99):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == "HanksBakun2002":
        try:
            # only small magnitude case
            assert a <= 537
        except AssertionError as e:
            e.args += (
                "Cannot use HanksAndBakun2002 equation for a > 537 km^2 (Which represents a rupture with 6.31 Mw)",)
            raise
        mw = np.log10(a) + 3.98

    elif mw_scaling_rel == "BerrymanEtAl2002":
        mw = np.log10(a) + 4.18
        # i.e. set L=W=sqrt(A) in Eqn 1 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    elif mw_scaling_rel == "VillamorEtAl2001":
        mw = np.log10(a) * 4 / 3 + 3.39
        # i.e. Eqn 2 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
        # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.

    elif mw_scaling_rel == "Leonard2014":
        if round(rake % 360 / 90.0) % 2:
            mw = np.log10(a) + leonard_ds
        else:
            mw = np.log10(a) + leonard_ss

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # magnitude
    return mw


def mag2mom(mw):
    """Converts magnitude to moment"""
    return np.exp(1.5 * (mw + 10.7) * np.log(10.0))


def mom2mag(mom):
    """Converts moment to magnitude"""
    return (2 / 3.0 * np.log(mom) / np.log(10.0)) - 10.7

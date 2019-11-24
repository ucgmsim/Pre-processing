from enum import Enum
from typing import Union

import numpy as np


class MagnitudeScalingRelations(Enum):
    HANKSBAKUN2002 = "hanksbakun2002"
    BERRYMANETAL2002 = "berrymanetal2002"
    VILLAMORETAL2001 = "villamoretal2001"
    LEONARD2014 = "leonard2014"


def mw_2_a_scaling_relation(
    mw: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[None, float] = None,
    leonard_ds: float = 4.00,
    leonard_ss: float = 3.99,
):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        A = mw_to_a_hanksbakum(mw)

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002:
        A = mw_to_a_berrymanetal(mw)

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        A = mw_to_a_villamoretal(mw)

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        A = mw_to_a_leonard(leonard_ds, leonard_ss, mw, rake)

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # Area
    return float(A)


def mw_to_a_hanksbakum(mw):
    if mw > 6.71:
        raise ValueError("Cannot use HanksAndBakun2002 equation for mw > 6.71")
    return 10 ** (mw - 3.98)


def mw_to_a_berrymanetal(mw):
    # i.e. set L=W=sqrt(A) in Eqn 1 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
    # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.
    return 10 ** (mw - 4.18)


def mw_to_a_villamoretal(mw):
    # i.e. Eqn 2 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
    # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.
    return 10 ** (0.75 * (mw - 3.39))


def mw_to_a_leonard(leonard_ds, leonard_ss, mw, rake):
    if round(rake % 360 / 90.0) % 2:
        A = 10 ** (mw - leonard_ds)
    else:
        A = 10 ** (mw - leonard_ss)
    return A


def a_2_mw_scaling_relation(
    a: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[float, None] = None,
    leonard_ds: float = 4.00,
    leonard_ss: float = 3.99,
):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        mw = a_to_mw_hanksbakum(a)

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002:
        mw = a_to_mw_berrymanetal(a)

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        mw = a_to_mw_villamoretal(a)

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        mw = a_to_mw_leonard(a, leonard_ds, leonard_ss, rake)

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # magnitude
    return float(mw)


def a_to_mw_leonard(a, leonard_ds, leonard_ss, rake):
    if round(rake % 360 / 90.0) % 2:
        mw = np.log10(a) + leonard_ds
    else:
        mw = np.log10(a) + leonard_ss
    return mw


def a_to_mw_villamoretal(a):
    # i.e. Eqn 2 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
    # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.
    return np.log10(a) * 4 / 3 + 3.39


def a_to_mw_berrymanetal(a):
    # i.e. set L=W=sqrt(A) in Eqn 1 in [1] Stirling, MW, Gerstenberger, M, Litchfield, N, McVerry, GH, Smith, WD, Pettinga, JR, Barnes, P. 2007.
    # updated probabilistic seismic hazard assessment for the Canterbury region, GNS Science Consultancy Report 2007/232, ECan Report Number U06/6. 58pp.
    return np.log10(a) + 4.18


def a_to_mw_hanksbakum(a):
    if a > 10 ** (6.71 - 3.98):
        raise ValueError(
            "Cannot use HanksAndBakun2002 equation for a > 537 km^2 (Which represents a rupture with 6.71 Mw)"
        )
    return np.log10(a) + 3.98


def mag2mom(mw):
    """Converts magnitude to moment"""
    return np.exp(3.0 / 2.0 * (mw + 10.7) * np.log(10.0))


def mom2mag(mom):
    """Converts moment to magnitude"""
    return (2.0 / 3.0 * np.log(mom) / np.log(10.0)) - 10.7

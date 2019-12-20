from enum import Enum
from typing import Union

import numpy as np


class MagnitudeScalingRelations(Enum):
    HANKSBAKUN2002 = "HANKSBAKUN2002"
    BERRYMANETAL2002 = "BERRYMANETAL2002"
    VILLAMORETAL2001 = "VILLAMORETAL2001"
    LEONARD2014 = "LEONARD2014"


def mw_2_lw_scaling_relation(
    mw: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[None, float] = None,
):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        l = w = np.sqrt(mw_to_a_hanksbakun(mw))

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002:
        l = w = np.sqrt(mw_to_a_berrymanetal(mw))

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        l = w = np.sqrt(mw_to_a_villamoretal(mw))

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        l, w = mw_to_lw_leonard(mw, rake)

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # Area
    return l, w


def mw_to_a_hanksbakun(mw):
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


def mw_to_lw_leonard(mw, rake):
    if round(rake % 360 / 90.0) % 2:
        A = 10 ** (mw - 4.0)
        l = 10 ** ((mw - 4.0) / 2)
        if l > 5.4:
            l = 10 ** ((mw - 4.24) / 1.667)
        w = 10 ** ((mw - 3.63) / 2.5)

    else:
        A = 10 ** (mw - 3.99)
        l = 10 ** ((mw - 4.17) / 1.667)
        if l > 45:
            l = 10 ** (mw - 5.27)
        w = 10 ** ((mw - 3.88) / 2.5)

    r = max(l / w, 1)
    w = np.sqrt(A / r)
    l = r * w

    return l, w


def wl_to_mw_leonard(l, w, rake):
    return a_to_mw_leonard(w * l, 4.00, 3.99, rake)


def lw_2_mw_scaling_relation(
    l: float,
    w: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[float, None] = None,
):
    """
    Return the fault Area from the mw and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        mw = a_to_mw_hanksbakun(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002:
        mw = a_to_mw_berrymanetal(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        mw = a_to_mw_villamoretal(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        mw = wl_to_mw_leonard(l, w, rake)

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


def a_to_mw_hanksbakun(a):
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

from enum import Enum
from typing import Union

import numpy as np


MAGNITUDE_ROUNDING_THRESHOLD = 7.5


class MagnitudeScalingRelations(Enum):
    HANKSBAKUN2002 = "HANKSBAKUN2002"
    BERRYMANETAL2002 = "BERRYMANETAL2002"
    VILLAMORETAL2001 = "VILLAMORETAL2001"
    LEONARD2014 = "LEONARD2014"
    SKARLATOUDIS2016 = "SKARLATOUDIS2016"
    STIRLING2008 = "STIRLING2008"


def get_area(fault: "Fault"):
    if fault.magnitude_scaling_relation == MagnitudeScalingRelations.HANKSBAKUN2002:
        farea = mw_to_a_hanksbakun(fault.magnitude)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.BERRYMANETAL2002:
        farea = mw_to_a_berrymanetal(fault.magnitude)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.VILLAMORETAL2001:
        farea = mw_to_a_villamoretal(fault.magnitude)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.LEONARD2014:
        farea = mw_to_a_leonard(fault.magnitude, fault.rake)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.SKARLATOUDIS2016:
        farea = mw_to_a_skarlatoudis(fault.magnitude)

    else:
        raise ValueError(
            "Invalid mw_scaling_rel: {}. Exiting.".format(
                fault.magnitude_scaling_relation
            )
        )

    # Area
    return farea


def get_width(fault: "Fault"):
    if fault.magnitude_scaling_relation == MagnitudeScalingRelations.HANKSBAKUN2002:
        fwidth = np.sqrt(mw_to_a_hanksbakun(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.BERRYMANETAL2002:
        fwidth = np.sqrt(mw_to_a_berrymanetal(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.VILLAMORETAL2001:
        fwidth = np.sqrt(mw_to_a_villamoretal(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.LEONARD2014:
        fwidth = mw_to_w_leonard(fault.magnitude, fault.rake)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.SKARLATOUDIS2016:
        fwidth = np.sqrt(mw_to_a_skarlatoudis(fault.magnitude))

    else:
        raise ValueError(
            "Invalid mw_scaling_rel: {}. Exiting.".format(
                fault.magnitude_scaling_relation
            )
        )
    return fwidth


def get_length(fault: "Fault"):
    if fault.magnitude_scaling_relation == MagnitudeScalingRelations.HANKSBAKUN2002:
        flength = np.sqrt(mw_to_a_hanksbakun(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.BERRYMANETAL2002:
        flength = np.sqrt(mw_to_a_berrymanetal(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.VILLAMORETAL2001:
        flength = np.sqrt(mw_to_a_villamoretal(fault.magnitude))

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.LEONARD2014:
        flength = mw_to_l_leonard(fault.magnitude, fault.rake)

    elif fault.magnitude_scaling_relation == MagnitudeScalingRelations.SKARLATOUDIS2016:
        flength = np.sqrt(mw_to_a_skarlatoudis(fault.magnitude))

    else:
        raise ValueError(
            "Invalid mw_scaling_rel: {}. Exiting.".format(
                fault.magnitude_scaling_relation
            )
        )
    return flength


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

    elif mw_scaling_rel == MagnitudeScalingRelations.SKARLATOUDIS2016:
        l = w = np.sqrt(mw_to_a_skarlatoudis(mw))

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


def mw_to_lw_leonard(mw, rake):
    """
    Calculates fault width and length from Leonards 2014 magnitude scaling relations
    :param mw: The magnitude of the fault or event
    :param rake: The rake of the event, for classifying strike slip or dip slip
    :return: The fault length and width
    """
    farea = mw_to_a_leonard(mw, rake)
    flength = mw_to_l_leonard(mw, rake)
    fwidth = mw_to_w_leonard(mw, rake)

    r = max(flength / fwidth, 1)
    fwidth = np.sqrt(farea / r)
    flength = r * fwidth

    return flength, fwidth


def mw_to_a_leonard(mw, rake):
    if round(rake % 360 / 90.0) % 2:
        farea = 10 ** (mw - 4.0)
    else:
        farea = 10 ** (mw - 3.99)
    return farea


def mw_to_l_leonard(mw, rake):
    if round(rake % 360 / 90.0) % 2:
        flength = 10 ** ((mw - 4.0) / 2)
        if flength > 5.4:
            flength = 10 ** ((mw - 4.24) / 1.667)
    else:
        flength = 10 ** ((mw - 4.17) / 1.667)
        if flength > 45:
            flength = 10 ** (mw - 5.27)
    return flength


def mw_to_w_leonard(mw, rake):
    if round(rake % 360 / 90.0) % 2:
        fwidth = 10 ** ((mw - 3.63) / 2.5)
    else:
        fwidth = 10 ** ((mw - 3.88) / 2.5)
    return fwidth


def wl_to_mw_leonard(l, w, rake):
    return a_to_mw_leonard(w * l, 4.00, 3.99, rake)


def lw_2_mw_scaling_relation(
    l: float,
    w: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[float, None] = None,
):
    """
    Return the mw from the fault Area and a mw Scaling relation.
    """
    if mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        mw = a_to_mw_hanksbakun(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.BERRYMANETAL2002:
        mw = a_to_mw_berrymanetal(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        mw = a_to_mw_villamoretal(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        mw = wl_to_mw_leonard(l, w, rake)

    elif mw_scaling_rel == MagnitudeScalingRelations.SKARLATOUDIS2016:
        mw = a_to_mw_skarlatoudis(l * w)

    elif mw_scaling_rel == MagnitudeScalingRelations.STIRLING2008:
        mw = lw_2_mw_stirling(l, w)

    else:
        raise ValueError("Invalid mw_scaling_rel: {}. Exiting.".format(mw_scaling_rel))

    # magnitude
    return float(mw)


def lw_2_mw_sigma_scaling_relation(
    l: float,
    w: float,
    mw_scaling_rel: MagnitudeScalingRelations,
    rake: Union[float, None] = None,
):
    """
    Return the mw and corresponding sigma from fault Area for a mw Scaling relation.
    """
    sigma = None
    # magnitude
    mw = lw_2_mw_scaling_relation(
        l,
        w,
        mw_scaling_rel,
        rake,
    )
    if mw_scaling_rel == MagnitudeScalingRelations.LEONARD2014:
        sigma = mw_sigma_leonard(rake)
    elif mw_scaling_rel == MagnitudeScalingRelations.SKARLATOUDIS2016:
        sigma = mw_sigma_skarlatoudis()
    elif mw_scaling_rel == MagnitudeScalingRelations.HANKSBAKUN2002:
        sigma = mw_sigma_hanksbakun()
    elif mw_scaling_rel == MagnitudeScalingRelations.VILLAMORETAL2001:
        sigma = mw_sigma_villamor()
    elif mw_scaling_rel == MagnitudeScalingRelations.STIRLING2008:
        sigma = mw_sigma_stirling()
    return mw, sigma


def mw_sigma_stirling():
    return 0.18


def mw_sigma_skarlatoudis():
    return np.log10(1.498)  # 1.498 is the sigma value from Table 3 in Skarlatoudis 2016


def mw_sigma_leonard(rake):
    """
    sigma values are determined from the Leonard 2014 paper, based on Table 3 - S(a) passing it through the equation

    rake: to determine if DS or SS
    returns sigma value for mean Mw
    """
    if round(rake % 360 / 90.0) % 2:
        sigma = 0.3
    else:
        sigma = 0.26
    return sigma


def mw_sigma_hanksbakun():
    return 0.22


def mw_sigma_villamor():
    return 0.195


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
        return 3.07 + (4.0 / 3.0) * np.log(a) * 0.434294
    return np.log10(a) + 3.98


def a_to_mw_skarlatoudis(a):
    return np.log10(a) - (np.log10(1.77 * np.power(10.0, -10)) + 6.03)


def mw_to_a_skarlatoudis(mw):
    return 1.77 * 10 ** -10 * (10 ** ((3 * (mw + 6.03)) / 2)) ** (2 / 3)


def lw_2_mw_stirling(l, w):
    return 4.18 + (2.0 / 3.0) * np.log10(w) + (4.0 / 3.0) * np.log10(l)


def mag2mom(mw):
    """Converts magnitude to moment - dyne-cm"""
    return 10 ** (3.0 / 2.0 * (mw + 10.7))


def mom2mag(mom):
    """Converts moment to magnitude - dyne-cm"""
    return (2.0 / 3.0 * np.log10(mom)) - 10.7


def mag2mom_nm(mw):
    """Converts magnitude to moment - newtonmetre"""
    return 10 ** (9.05 + 1.5 * mw)


def mom2mag_nm(mom):
    """Converts moment to magnitude - newtonmetre"""
    return (np.log10(mom) - 9.05) / 1.5


def round_subfault_size(dist, mag):
    if mag > MAGNITUDE_ROUNDING_THRESHOLD:
        return round(dist * 2) / 2
    else:
        return round(dist * 10) / 10

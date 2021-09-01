import numpy as np


def mw_to_a_strasser_2010_interface(mw):
    """
    strasser2010scaling
    interface [6.3-9.4]
    """

    if mw < 6.3 or mw > 9.4:
        raise ValueError("Magnitude out of range for model")

    area = np.power(10, 0.952 * mw - 3.476)
    return area


def a_to_mw_strasser_2010_interface(a):
    """
    Inversion of the above equation
    strasser2010scaling
    interface [6.3-9.4]
    """

    if a < 332 or a > 297000:
        raise ValueError("Magnitude out of range for model")

    mw = (np.log10(a) + 3.476) / 0.952
    return mw


def mw_to_w_strasser_2010_interface(mw):
    if mw < 6.3 or mw > 9.4:
        raise ValueError("Magnitude out of range for model")

    w = np.power(10, 0.351 * mw - 0.882)
    return w


def mw_to_l_strasser_2010_interface(mw):
    if mw < 6.3 or mw > 9.4:
        raise ValueError("Magnitude out of range for model")

    l = np.power(10, 0.585 * mw - 2.477)
    return l


def mw_to_a_strasser_2010_slab(mw):
    """
    strasser2010scaling
    slab [5.9-7.8]
    """

    if mw < 5.9 or mw > 7.8:
        raise ValueError("Magnitude out of range for model")

    area = np.power(10, 0.890 * mw - 3.225)
    return area


def mw_to_w_strasser_2010_slab(mw):
    if mw < 5.9 or mw > 7.8:
        raise ValueError("Magnitude out of range for model")

    w = np.power(10, 0.356 * mw - 1.058)
    return w


def mw_to_l_strasser_2010_slab(mw):
    if mw < 5.9 or mw > 7.8:
        raise ValueError("Magnitude out of range for model")

    l = np.power(10, 0.562 * mw - 2.350)
    return l

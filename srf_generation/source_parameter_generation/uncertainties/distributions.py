import random

import numpy as np
from scipy.stats import truncnorm


def relative_uniform(mean, scale_factor):
    """Returns a value from a uniform distribution where the min and max
    are scaled relative to the middle value given as the mean"""
    return random.uniform(mean * (1 - scale_factor), mean * (1 + scale_factor))


def uniform(mean, half_range):
    return random.uniform(mean - half_range, mean + half_range)


def truncated_normal(mean, std_dev, std_dev_limit=2):
    return float(
        truncnorm(-std_dev_limit, std_dev_limit, loc=mean, scale=std_dev).rvs()
    )


def weibull(k=3.353, scale_factor=0.612):
    """Weibull distribution. Defaults are for nhm2srf dhypo generation"""
    return scale_factor * np.random.weibull(k)


def truncated_log_normal(mean, std_dev, std_dev_limit=2) -> float:
    return np.exp(
        truncnorm(
            -std_dev_limit, std_dev_limit, loc=np.log(float(mean)), scale=std_dev
        ).rvs()
    )


def rand_shyp():
    """Generates a value for the hypocentre along the length of the fault. Uses defaults from nhm2srf"""
    # normal distribution
    shyp_mu = 0.0
    shyp_sigma = 0.25

    shyp = truncated_normal(shyp_mu, shyp_sigma)
    return shyp

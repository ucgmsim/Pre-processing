from random import uniform

import numpy as np
from scipy.stats import truncnorm


def relative_uniform_distribution(mean, scale_factor):
    """Returns a value from a uniform distribution where the min and max
    are scaled relative to the middle value given as the mean"""
    return uniform(mean * (1 - scale_factor), mean * (1 + scale_factor))


def uniform_distribution(mean, half_range):
    return uniform(mean - half_range, mean + half_range)


def truncated_normal(mean, std_dev, std_dev_limit=2):
    return float(
        truncnorm(-std_dev_limit, std_dev_limit, loc=mean, scale=std_dev).rvs()
    )


def weibull(k=3.353, scale_factor=0.612):
    return scale_factor * np.random.weibull(k)


def truncated_log_normal(mean, std_dev, std_dev_limit=2) -> float:
    return np.exp(
        truncnorm(
            -std_dev_limit, std_dev_limit, loc=np.log(float(mean)), scale=std_dev
        ).rvs()
    )

from random import uniform, normalvariate

import numpy as np
from scipy.stats import truncnorm


def relative_uniform_distribution(mean, scale_factor):
    """A uniform distribution where the min and max are scaled relative to the middle value given as the mean"""
    return uniform(mean * (1 - scale_factor), mean * (1 + scale_factor))


def uniform_distribution(mean, half_range):
    return uniform(mean - half_range, mean + half_range)


def truncated_normal(mean, std_dev, std_dev_limit=2):
    return float(
        truncnorm(-std_dev_limit, std_dev_limit, loc=mean, scale=std_dev).rvs()
    )


def weibull(k=3.353, scale_factor=0.612):
    return scale_factor * np.random.weibull(k)


def duplicate_value(targets, *perturbated_dicts):
    for p_dict in perturbated_dicts:
        if targets in p_dict.keys():
            return p_dict[targets]
    raise ValueError("Unable to find target to duplicate")


def trunc_log_normal(mean, std_dev, std_dev_limit=2):
    return np.exp(
        truncnorm(
            -std_dev_limit, std_dev_limit, loc=np.log(float(mean[4:-1])), scale=std_dev
        ).rvs()
    )

#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[14]:


"""This version includes the following revisions from v800: 

- vs profile perturbation removed
- changed Qs LF and HF perturbation to be correlated
- changed Qs std_dev_{ln} to be 0.5 (previously 0.3)
- updated R-distance and depth uncertainties to be event specific
- using a normal distribution for depth uncertainty
- removed kappa uncertainty
- added site model uncertainty to Vs30 uncertainty

"""

from typing import Any, Dict, Union

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv
import math
import matplotlib.lines as mlines
from numpy.random import randn
from pathlib import Path

import os
from scipy.stats import norm
from scipy.stats import rv_continuous
from scipy.stats import kstest
import matplotlib.style
import matplotlib as mpl
from time import time
from pandas import DataFrame
import qcore
from qcore.uncertainties import distributions
from qcore import geo
from typing import Any, Dict, Union  # new line
from geographiclib.geodesic import Geodesic

from srf_generation.Fault import fault_factory, Type1

from srf_generation.source_parameter_generation.common import get_depth_property
from srf_generation.source_parameter_generation.uncertainties.common import (
    verify_realisation_params,
    GCMT_Source,
)
from qcore.uncertainties.mag_scaling import (
    lw_to_mw_sigma_scaling_relation,
    MagnitudeScalingRelations,
    mag2mom,
)


def verify_1d_vel_mod(vel_mod_1d: pd.DataFrame):  # new line
    for line in vel_mod_1d.itertuples():  # new line
        pass
        # assert line._fields == ('depth', 'vp', 'vs', 'rho', 'qp', 'qs')       #new line


def generate_source_params(
    source_data: Union[GCMT_Source],
    additional_source_parameters: Dict[str, Any],  # new line
    vel_mod_1d: pd.DataFrame = None,  # new line
    vs30_data: pd.DataFrame = None,
    r_rup_data: pd.DataFrame = None,
    **kwargs,  # new line
) -> Dict[str, Any]:
    """source_data should have the following parameters available via . notation:
    - source_data.pid: name of the event
    - source_data.lat: latitude
    - source_data.lon: longitude
    - source_data.depth: depth
    - source_data.mag: magnitude
    - source_data.strike
    - source_data.dip
    - source_data.rake
    """

    print(source_data.pid)

    realisation = kwargs  # new line
    ### Start of custom code area

    qs_z = distributions.truncated_normal(0, 1, 2.5)

    if isinstance(source_data, GCMT_Source):
        realisation["params"] = generate_from_gcmt(
            source_data, qs_z, additional_source_parameters
        )  # modified
    else:
        raise ValueError("Not a valid type")

    if vel_mod_1d is not None:  # new line
        realisation["vel_mod_1d"] = vel_mod_1d  # new line

    # onedprofile_ori = pd.read_csv("/nesi/project/nesi00213/Pre-processing/SrfGen/lp_generic1d-gp01_v1.vmod",delim_whitespace=True, names=['depth', 'vp', 'vs', 'rho', 'qp', 'qs'],skiprows=1)
    # onedprofile_pert = onedprofile_ori.copy()
    # vs_profile = [distributions.truncated_log_normal(vs, 0.05, 4) for vs in onedprofile_pert.vs]
    # onedprofile_pert.vs = vs_profile

    # realisation["vel_mod_1d"] = onedprofile_pert
    # if realisation["params"]["type"] == 1:
    # realisation["params"]["vs"] = get_depth_property(source_data.depth, onedprofile_pert, "vs")

    # HF 1d profile
    HF_onedprofile_ori = pd.read_csv(
        "/nesi/project/nesi00213/VelocityModel/Mod-1D/Cant1D_v3-midQ_OneRay.1d",
        delim_whitespace=True,
        names=["depth", "vp", "vs", "rho", "qp", "qs"],
        skiprows=1,
    )
    HF_onedprofile_pert = HF_onedprofile_ori.copy()

    # qs_pert = distributions.truncated_log_normal(1, 0.3, 2.5)         #********
    # HF_qs_profile = HF_onedprofile_pert.qs*qs_pert                    #********

    qs_pert = np.exp(np.log(1) + qs_z * 0.5)

    # HF_qs_profile = [distributions.truncated_log_normal(qs, 0.3, 2.5) for qs in HF_onedprofile_pert.qs]
    #HF_qs_profile = [qs for qs in HF_onedprofile_pert.qs]

    # qs_z = (HF_qs_profile - HF_onedprofile_pert.qs)/0.3    #*****

    #HF_onedprofile_pert.qs = HF_qs_profile

    realisation["hf_vel_mod_1d"] = HF_onedprofile_pert

    #measurement_uncertainty = 0.1
    #site_model_uncertainty = 0.3
    #total_vs30_uncertainty = vs30_data["sigma"].values
    total_vs30_uncertainty = 0.1

    if vs30_data is not None:
        print("Got vs30")
        realisation["vs30"] = vs30_data.copy(deep=True)
        #print(vs30_data["median"].values)
        realisation["vs30"]["vs30"] = distributions.truncated_log_normal(vs30_data["median"].values, total_vs30_uncertainty, 4)
    else:
        print("Didn't get vs30")

    #if realisation["vs30"]["vs30"] <= 0:
    #    realisation["vs30"]["vs30"] = 0






    ### End of custom code area

    # verify_realisation_params(realisation["params"])                      #modified
    # realisation['params']['qs_z'] = ' '.join(list(qs_z.astype(str)))    #*****
    realisation["params"]["qs_pert"] = qs_pert  # ********

    if vel_mod_1d is not None:  # new
        verify_1d_vel_mod(realisation["vel_mod_1d"])  # new

    return realisation  # modified


# distribution functions:


def uniform_dist(u_mean, u_half_range):
    """function for sampling uniform distribution"""
    u_sample = np.random.uniform(
        low=u_mean - u_half_range, high=u_mean + u_half_range, size=1
    )
    return u_sample


def generate_from_gcmt(
    sources_line: GCMT_Source, qs_z, additional_source_parameters: Dict[str, Any]
):  # modified
    # area = mw_2_a_scaling_relation(
    #    sources_line.mag,
    #    MagnitudeScalingRelations.LEONARD2014.value,
    #    sources_line.strike,
    # )

    #rdist_sigma = 2
    #depth_sigma = additional_source_parameters["depth_sigma"]

    sites_info_csv = Path(additional_source_parameters["sites_info"])
    assert sites_info_csv.exists()
    sites_info = pd.read_csv(sites_info_csv)


    ### the following parameters feed into srfgen

    mag = distributions.truncated_normal(sources_line.mag, 0.2, 4) # magnitude
    # rvfrac = uniform_dist(0.8, 0.075)                   #rupture velocity factor

    # lat_temp,lon_temp = geo.ll_shift(sources_line.lat, sources_line.lon, distributions.truncated_normal(0.0, 1.0, 2), 0)
    # lat,lon = geo.ll_shift(lat_temp, lon_temp, distributions.truncated_normal(0.0, 1.0, 2), 90)

    #theta = uniform_dist(0, 180)
    # r_distance = (distributions.truncated_normal(0.0, 3.65, 2))*1000
    #r_distance = (distributions.truncated_normal(0.0, rdist_sigma, 4)) * 1000

    if sites_info is not None:
        rrup_data = sites_info.rrup
        rrup_data_perturbed = distributions.truncated_log_normal(rrup_data , 0.1, 4)
        sites_info["rrup"] = rrup_data_perturbed
        new_sites_info_csv = sites_info_csv.parent / f"{sites_info_csv.stem}_rrup_perturbed.csv"
        sites_info.to_csv(new_sites_info_csv, index=False)
    else:
        print("Didn't get sites_info CSV")


    #temp = Geodesic.WGS84.Direct(sources_line.lat, sources_line.lon, theta, r_distance)
    #lat = distributions.truncated_normal(sources_line.lat, 1, 4)
    lat = sources_line.lat
    #lon = distributions.truncated_normal(sources_line.lon, 1, 4)
    lon = sources_line.lon

    # depth = distributions.truncated_log_normal(sources_line.depth, 0.3, 2)

    #depth = distributions.truncated_normal(sources_line.depth, depth_sigma, 2)
    depth = sources_line.depth
    #if depth <= 3:
    #    depth = 3
    #if depth >= 37:
    #    depth = 37

    #strike = distributions.truncated_normal(sources_line.strike, 10, 2)
    strike = sources_line.strike
    #dip = distributions.truncated_normal(sources_line.dip, 10, 2)
    dip = sources_line.dip
    #rake = distributions.truncated_normal(sources_line.rake, 15, 4)
    rake = sources_line.rake

    # qsfrac = distributions.truncated_log_normal(50, 0.3, 2.5)
    qsfrac = 50

    #unmod_strike = strike
    #unmod_dip = dip
    #unmod_rake = rake

    # calculate normalised perturbations:
    #z_mw_pert = (mag - sources_line.mag) / 0.075
    #z_lat_pert = (
    #    geo.ll_dist(sources_line.lon, sources_line.lat, sources_line.lon, lat) / 1.0
    #)
    #z_long_pert = (
    #    geo.ll_dist(sources_line.lon, sources_line.lat, lon, sources_line.lat) / 1.0
    #)
    #z_depth_pert = (depth - sources_line.depth) / 2.0
    #z_strike_pert = (strike - sources_line.strike) / 10.0
    #z_dip_pert = (dip - sources_line.dip) / 10.0
    #z_rake_pert = (rake - sources_line.rake) / 15.0

    # correct strike dip and rake:
    #if dip > 90:
    #    strike = (strike + 180) % 360  # flip the strike
    #    dip = 180 - dip
    #    rake = rake * -1

    #if dip < 0:
    #    strike = (strike + 180) % 360
    #    dip = 90 - dip % 90
    #    rake = rake * -1

    #if rake > 180:
    #    rake = rake - 360

    #if rake < -180:
    #    rake = rake + 360

    #if strike > 360:
    #    strike = strike % 360

    #if strike < 0:
    #    strike = strike % 360

    ### the following parameters feed into the sim_params.yaml

    sdrop = 50
    rvfac = 0.8  # rupture velocity factor
    if rvfac == 0:
        rvfac = 0.001
    # kappa = distributions.truncated_log_normal(0.045, 0.3, 2)
    kappa = 0.045

    dpath_pert = 0

    # calculate normalised perturbations:
    #z_sdrop_pert = (sdrop - 50) / 1.35
    #z_rvfac_pert = (rvfac - 0.8) / 0.075
    # z_kappa_pert = (kappa - 0.045)/1.35

    params = {
        "type": 1,
        "name": sources_line.pid,
        "latitude": lat,
        "longitude": lon,
        "depth": depth,
        "magnitude": mag,
        "moment": float(mag2mom(mag)),
        "strike": strike,
        "dip": dip,
        "rake": rake,
        "kappa": kappa,
        "sdrop": sdrop,
        "rvfac": rvfac,
        "qsfrac": qsfrac,
        "dpath_pert": dpath_pert
    }

    params.update(additional_source_parameters)  # modified
    return params

import argparse
from subprocess import call, PIPE
from typing import Dict, Any, Union
from tempfile import NamedTemporaryFile

import yaml
from h5py import File as h5open
import numpy as np
from os import makedirs, path, remove
import pandas as pd
from qcore import binary_version, srf, geo

from srf_generation.source_parameter_generation.uncertainties.common import (
    HF_RUN_PARAMS,
    BB_RUN_PARAMS,
    LF_RUN_PARAMS,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import mag2mom

SRF2STOCH = "srf2stoch"
GENERICSLIP2SRF = "generic_slip2srf"
FAULTSEG2GSFDIPDIR = "fault_seg2gsf_dipdir"


def gen_stoch(stoch_file, srf_file, single_segment=False):
    """
    Creates stoch file from srf file.
    """
    out_dir = path.dirname(stoch_file)
    makedirs(out_dir, exist_ok=True)
    dx, dy = 2.0, 2.0
    if not srf.is_ff(srf_file):
        dx, dy = srf.srf_dxy(srf_file)
    with open(stoch_file, "w") as stochp, open(srf_file, "r") as srfp:
        if single_segment:
            call(
                [
                    binary_version.get_unversioned_bin(SRF2STOCH),
                    f"target_dx={dx}",
                    f"target_dy={dy}",
                ],
                stdin=srfp,
                stdout=stochp,
            )
        else:
            call(
                [binary_version.get_unversioned_bin(SRF2STOCH), f"dx={dx}", f"dy={dy}"],
                stdin=srfp,
                stdout=stochp,
            )


def create_info_file(
    srf_file,
    srf_type,
    mag,
    rake,
    dt,
    vm=None,
    vs=None,
    rho=None,
    centroid_depth=None,
    lon=None,
    lat=None,
    shypo=None,
    dhypo=None,
    mwsr=None,
    tect_type=None,
    dip_dir=None,
    file_name=None,
):
    """
    Stores SRF metadata as hdf5.
    srf_file: SRF path used as basename for info file and additional metadata
    """
    planes = srf.read_header(srf_file, idx=True)
    hlon, hlat, hdepth = srf.get_hypo(srf_file, depth=True)

    dbottom = []
    corners = np.zeros((len(planes), 4, 2))
    for i, p in enumerate(planes):
        # currently only support single dip dir value
        if dip_dir is not None:
            dip_deg = dip_dir
        else:
            dip_deg = p["strike"] + 90
        # projected fault width (along dip direction)
        pwid = p["width"] * np.cos(np.radians(p["dip"]))
        corners[i, 0] = geo.ll_shift(
            p["centre"][1], p["centre"][0], p["length"] / 2.0, p["strike"] + 180
        )[::-1]
        corners[i, 1] = geo.ll_shift(
            p["centre"][1], p["centre"][0], p["length"] / 2.0, p["strike"]
        )[::-1]
        corners[i, 2] = geo.ll_shift(corners[i, 1, 1], corners[i, 1, 0], pwid, dip_deg)[
            ::-1
        ]
        corners[i, 3] = geo.ll_shift(corners[i, 0, 1], corners[i, 0, 0], pwid, dip_deg)[
            ::-1
        ]
        dbottom.append(p["dtop"] + p["width"] * np.sin(np.radians(p["dip"])))

    if file_name is None:
        file_name = srf_file.replace(".srf", ".info")

    with h5open(file_name, "w") as h:
        a = h.attrs
        # only taken from given parameters
        a["type"] = srf_type
        a["dt"] = dt
        a["rake"] = rake
        a["mag"] = mag
        # srf header data
        for k in planes[0].keys():
            a[k] = [p[k] for p in planes]
        # point source has 1 vs/rho, others are taken from vm
        if srf_type == 1:
            a["vs"] = vs
            a["rho"] = rho
            a["hlon"] = lon
            a["hlat"] = lat
            a["hdepth"] = centroid_depth
        else:
            a["vm"] = np.string_(path.basename(vm))
            a["hlon"] = hlon
            a["hlat"] = hlat
            a["hdepth"] = hdepth
            a["shyp0"] = shypo
            a["dhyp0"] = dhypo
        if mwsr is not None:
            a["mwsr"] = np.string_(mwsr)
        if tect_type is not None:
            a["tect_type"] = tect_type
        # derived parameters
        a["corners"] = corners
        a["dbottom"] = dbottom


def create_ps_srf(
    realisation_file: str,
    parameter_dictionary: Dict[str, Any],
    stoch_file: Union[None, str] = None,
):
    latitude = parameter_dictionary.pop("latitude")
    longitude = parameter_dictionary.pop("longitude")
    depth = parameter_dictionary.pop("depth")
    magnitude = parameter_dictionary.pop("magnitude")
    strike = parameter_dictionary.pop("strike")
    rake = parameter_dictionary.pop("rake")
    dip = parameter_dictionary.pop("dip")
    moment = parameter_dictionary.pop("moment", 0)
    dt = parameter_dictionary.pop("dt", 0.005)
    vs = parameter_dictionary.pop("vs", 3.20)
    rho = parameter_dictionary.pop("rho", 2.44)
    target_area_km = parameter_dictionary.pop("target_area_km", None)
    target_slip_cm = parameter_dictionary.pop("target_slip_cm", None)
    stype = parameter_dictionary.pop("stype", "cos")
    risetime = parameter_dictionary.pop("risetime", 0.5)
    inittime = parameter_dictionary.pop("inittime", 0.0)

    if moment <= 0:
        moment = mag2mom(magnitude)

    # size (dd) and slip
    if target_area_km is not None:
        dd = np.sqrt(target_area_km)
        slip = (moment * 1.0e-20) / (target_area_km * vs * vs * rho)
    elif target_slip_cm is not None:
        dd = np.sqrt(moment * 1.0e-20) / (target_slip_cm * vs * vs * rho)
        slip = target_slip_cm
    else:
        aa = np.exp(2.0 / 3.0 * np.log(moment) - 14.7 * np.log(10.0))
        dd = np.sqrt(aa)
        slip = (moment * 1.0e-20) / (aa * vs * vs * rho)

    ###
    ### file names
    ###
    srf_file = realisation_file.replace(".csv", ".srf")
    gsf_file = NamedTemporaryFile(mode='w', delete=False)

    ###
    ### create GSF
    ###
    with gsf_file as gsfp:
        gsfp.write("# nstk= 1 ndip= 1\n")
        gsfp.write(f"# flen= {dd:10.4f} fwid= {dd:10.4f}\n")
        gsfp.write(
            "# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n"
        )
        gsfp.write("1\n")
        gsfp.write(
            f"{longitude:11.5f} {latitude:11.5f} {depth:8.4f} {dd:8.4f} {dd:8.4f} "
            f"{strike:6.1f} {dip:6.1f} {rake:6.1f} {slip:8.2} {inittime:8.3f}    0\n"
        )

    ###
    ### create SRF
    ###
    commands = [
        binary_version.get_unversioned_bin(GENERICSLIP2SRF),
        f"infile={gsf_file.name}",
        f"outfile={srf_file}",
        "outbin=0",
        f"stype={stype}",
        f"dt={dt}",
        "plane_header=1",
        f"risetime={risetime}",
        "risetimefac=1.0",
        "risetimedep=0.0",
    ]
    call(commands, stderr=PIPE)

    remove(gsf_file)

    ###
    ### save STOCH
    ###
    if stoch_file is None:
        stoch_file = realisation_file.replace(".csv", ".stoch")
    gen_stoch(stoch_file, srf_file)

    ###
    ### save INFO
    ###
    create_info_file(
        srf_file,
        1,
        magnitude,
        rake,
        dt,
        vs=vs,
        rho=rho,
        centroid_depth=depth,
        lon=longitude,
        lat=latitude,
    )


def generate_sim_params_yaml(sim_params_file: str, parameters: Dict[str, Any]):
    makedirs(path.dirname(sim_params_file), exist_ok=True)

    sim_params = {}
    for key, value in parameters.items():
        if key in HF_RUN_PARAMS:
            if "hf" not in sim_params.keys():
                sim_params["hf"] = {}
            sim_params["hf"][key] = value

        elif key in BB_RUN_PARAMS:
            if "bb" not in sim_params.keys():
                sim_params["bb"] = {}
            sim_params["bb"][key] = value

        elif key in LF_RUN_PARAMS:
            if "emod3d" not in sim_params.keys():
                sim_params["emod3d"] = {}
            sim_params["emod3d"][key] = value

        else:
            sim_params[key] = value

    with open(sim_params_file, "w") as spf:
        yaml.dump(sim_params, spf)


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("realisation_file", type=path.abspath)
    return parser.parse_args()


def main():
    args = load_args()
    rel_df: pd.DataFrame = pd.read_csv(args.realisation_file)
    realisation = rel_df.to_dict(orient="records")[0]
    if realisation["type"] == 1:
        create_ps_srf(args.realisation_file, realisation)
    else:
        raise ValueError(
            f"Type {realisation['type']} faults are not currently supported. "
            f"Contact the software team if you believe this is an error."
        )
    sim_params_file = args.realisation_file.replace(".csv", ".yaml")
    generate_sim_params_yaml(sim_params_file, realisation)


if __name__ == "__main__":
    main()

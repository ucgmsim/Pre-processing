import argparse
from logging import Logger
from subprocess import call, PIPE, Popen
from typing import Dict, Any, Union
from tempfile import NamedTemporaryFile

import yaml
from h5py import File as h5open
import numpy as np
from os import makedirs, path
import pandas as pd
from qcore import binary_version, srf, geo, qclogging
from qcore.utils import compare_versions

from srf_generation.pre_processing_common import calculate_corners, get_hypocentre
from srf_generation.source_parameter_generation.uncertainties.common import (
    HF_RUN_PARAMS,
    BB_RUN_PARAMS,
    LF_RUN_PARAMS,
    get_seed,
)
from srf_generation.source_parameter_generation.uncertainties.mag_scaling import mag2mom

SRF2STOCH = "srf2stoch"
GENERICSLIP2SRF = "generic_slip2srf"
FAULTSEG2GSFDIPDIR = "fault_seg2gsf_dipdir"

CORNERS_HEADER = (
    "> header line here for specifics \n\
> This is the standard input file format where the \
hypocenter is first then for each \n\
>Hypocenter (reference??) \n",
    "> Below are the corners \
(first point repeated as fifth to close box \n",
)


def create_stoch(
    stoch_file,
    srf_file,
    single_segment=False,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Creates stoch file from srf file.
    """
    logger.debug("Generating stoch file")
    out_dir = path.dirname(stoch_file)
    makedirs(out_dir, exist_ok=True)
    dx, dy = 2.0, 2.0
    if not srf.is_ff(srf_file):
        dx, dy = srf.srf_dxy(srf_file)
    logger.debug(f"Saving stoch to {stoch_file}")
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
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Stores SRF metadata as hdf5.
    srf_file: SRF path used as basename for info file and additional metadata
    """
    logger.debug("Generating srf info file")
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
    logger.debug("Saving info file to {}".format(file_name))
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
    logger: Logger = qclogging.get_basic_logger,
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

    logger.debug(
        "All srf generation parameters successfully obtained from the realisation file"
    )

    if moment <= 0:
        logger.debug("moment is negative, calculating from magnitude")
        moment = mag2mom(magnitude)

    # size (dd) and slip
    if target_area_km is not None:
        logger.debug(
            f"target_area_km given ({target_area_km}), using it to calculate fault edge length and slip"
        )
        dd = np.sqrt(target_area_km)
        slip = (moment * 1.0e-20) / (target_area_km * vs * vs * rho)
    elif target_slip_cm is not None:
        logger.debug(
            f"target_slip_cm given ({target_slip_cm}), using it to calculate fault edge length"
        )
        dd = np.sqrt(moment * 1.0e-20) / (target_slip_cm * vs * vs * rho)
        slip = target_slip_cm
    else:
        aa = np.exp(2.0 / 3.0 * np.log(moment) - 14.7 * np.log(10.0))
        dd = np.sqrt(aa)
        slip = (moment * 1.0e-20) / (aa * vs * vs * rho)
    logger.debug(f"Slip: {slip}, fault plane edge length {dd}")

    ###
    ### file names
    ###
    srf_file = realisation_file.replace(".csv", ".srf")
    logger.debug(f"Srf will be saved to {srf_file}")

    ###
    ### create GSF
    ###
    with NamedTemporaryFile(mode="w+") as gsfp:
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
        gsfp.flush()

        ###
        ### create SRF
        ###
        commands = [
            binary_version.get_unversioned_bin(GENERICSLIP2SRF),
            f"infile={gsfp.name}",
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

    ###
    ### save STOCH
    ###
    if stoch_file is None:
        stoch_file = realisation_file.replace(".csv", ".stoch")
    create_stoch(stoch_file, srf_file, single_segment=True)

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


def create_ps_ff_srf(
    realisation_file: str,
    parameter_dictionary: Dict[str, Any],
    stoch_file: Union[None, str] = None,
    logger: Logger = qclogging.get_basic_logger,
):
    # pops
    latitude = parameter_dictionary.pop("latitude")
    longitude = parameter_dictionary.pop("longitude")
    depth = parameter_dictionary.pop("depth")
    magnitude = parameter_dictionary.pop("magnitude")
    strike = parameter_dictionary.pop("strike")
    rake = parameter_dictionary.pop("rake")
    dip = parameter_dictionary.pop("dip")
    dt = parameter_dictionary.pop("dt", 0.005)
    flen = parameter_dictionary.pop("flen")
    dlen = parameter_dictionary.pop("dlen")
    fwid = parameter_dictionary.pop("fwid")
    dwid = parameter_dictionary.pop("dwid")
    dtop = parameter_dictionary.pop("dtop")
    shypo = parameter_dictionary.pop("shypo")
    dhypo = parameter_dictionary.pop("dhypo")
    mwsr = parameter_dictionary.pop("mwsr")
    genslip_version = parameter_dictionary.pop("genslip_version")
    slip_cov = parameter_dictionary.pop("slip_cov", 0.005)
    rough = parameter_dictionary.pop("rough", 0.005)
    seed = parameter_dictionary.pop("seed", get_seed())

    # gets
    vel_mod_1d = parameter_dictionary.get("vel_mod_1d")
    rvfac = parameter_dictionary.get("rvfac", 0.005)

    logger.debug(
        "All srf generation parameters successfully obtained from the realisation file"
    )

    gsf_file = NamedTemporaryFile(mode="w", delete=False)
    logger.debug(f"Gsf will be saved to the temporary file {gsf_file}")
    srf_file = realisation_file.replace(".csv", ".srf")
    logger.debug(f"Srf will be saved to {srf_file}")
    corners_file = realisation_file.replace(".csv", ".corners")
    logger.debug(f"The corners file will be saved to {corners_file}")

    logger.debug("Getting corners and hypocentre")
    corners = get_corners(latitude, longitude, flen, fwid, dip, strike)
    hypocentre = get_hypocentre(latitude, longitude, shypo, dhypo, strike, dip)
    logger.debug("Saving corners and hypocentre")
    write_corners(corners_file, hypocentre, corners)

    def get_n(fault_size, sub_fault_size):
        return f"{fault_size / sub_fault_size:0f}"

    nx = get_n(flen, dlen)
    ny = get_n(fwid, dwid)

    with NamedTemporaryFile(mode="w", delete=False) as gsfp:
        gen_gsf(
            gsfp,
            longitude,
            latitude,
            dtop,
            strike,
            dip,
            rake,
            flen,
            fwid,
            nx,
            ny,
            logger=logger,
        )
        gsfp.flush()
        gen_srf(
            srf_file,
            gsfp,
            magnitude,
            dt,
            nx,
            ny,
            seed,
            shypo,
            dhypo,
            vel_mod_1d,
            genslip_version=genslip_version,
            rvfac=rvfac,
            slip_cov=slip_cov,
            rough=rough,
            logger=logger,
        )

    if stoch_file is None:
        stoch_file = realisation_file.replace(".csv", ".stoch")
    create_stoch(stoch_file, srf_file, single_segment=True, logger=logger)

    # save INFO
    create_info_file(
        srf_file,
        2,
        magnitude,
        rake,
        dt,
        centroid_depth=depth,
        lon=None,
        lat=None,
        tect_type=None,
        dip_dir=None,
        mwsr=mwsr,
        shypo=shypo + 0.5 * flen,
        dhypo=dhypo,
        vm=vel_mod_1d,
        logger=logger,
    )


def get_corners(lat, lon, flen, fwid, dip, strike):
    """
    Return Corners of a fault. Indexes are as below.
    0 1
    2 3
    """
    lats, lons = calculate_corners(
        dip,
        np.array([-flen / 2.0, flen / 2.0]),
        np.array([-fwid, 0.0]),
        lat,
        lon,
        strike,
    )

    # joined (lon, lat) pairs
    return np.dstack((lons.flat, lats.flat))[0]


def write_corners(filename, hypocentre, corners):
    """
    Write a corners text file (used to plot faults).
    """
    with open(filename, "w") as cf:
        # lines beginning with '>' are ignored
        cf.write(CORNERS_HEADER[0])
        cf.write("{:f} {:f}\n".format(*hypocentre))
        cf.write(CORNERS_HEADER[1])
        # 0 1 - draw in order to close box
        # 2 3
        for i in [0, 1, 3, 2, 0]:
            cf.write("{:f} {:f}\n".format(*corners[i]))


def gen_srf(
    srf_file,
    gsf_file,
    magnitude,
    dt,
    nx,
    ny,
    seed,
    shypo,
    dhypo,
    velocity_model,
    genslip_version="3.3",
    rvfac=None,
    rough=None,
    slip_cov=None,
    logger: Logger = qclogging.get_basic_logger(),
):
    genslip_bin = binary_version.get_genslip_bin(genslip_version)
    if compare_versions(genslip_version, "5"):
        logger.debug(
            "Using genslip version {}. Using nx and ny".format(genslip_version)
        )
        xstk = "nx"
        ydip = "ny"
    else:
        logger.debug(
            "Using genslip version {}. Using nstk and ndip".format(genslip_version)
        )
        xstk = "nstk"
        ydip = "ndip"
    cmd = [
        genslip_bin,
        "read_erf=0",
        "write_srf=1",
        "read_gsf=1",
        "write_gsf=0",
        f"infile={gsf_file.name}",
        f"mag={magnitude}",
        f"{xstk}={nx}",
        f"{ydip}={ny}",
        "ns=1",
        "nh=1",
        f"seed={seed}",
        f"velfile={velocity_model}",
        f"shypo={shypo}",
        f"dhypo={dhypo}",
        f"dt={dt}",
        "plane_header=1",
        "srf_version=1.0",
    ]
    if rvfac is not None:
        cmd.append(f"rvfrac={rvfac}")
    if rough is not None:
        cmd.append(f"alpha_rough={rough}")
    if slip_cov is not None:
        cmd.append(f"slip_sigma={slip_cov}")
    logger.info("Creating SRF with command: {}".format(" ".join(cmd)))
    with open(srf_file, "w") as srfp:
        call(cmd, stdout=srfp)


def gen_gsf(
    gsfp,
    lon,
    lat,
    dtop,
    strike,
    dip,
    rake,
    flen,
    fwid,
    nx,
    ny,
    logger: Logger = qclogging.get_basic_logger(),
):
    logger.debug(f"Saving gsf to {gsfp.name}")
    gexec = Popen(
        [binary_version.get_unversioned_bin(FAULTSEG2GSFDIPDIR), "read_slip_vals=0"],
        stdin=PIPE,
        stdout=gsfp,
    )
    gexec.communicate(
        f"1\n{lon:f} {lat:f} {dtop:f} {strike} {dip} {rake} {flen:f} {fwid:f} {nx:s} {ny:s}".encode(
            "utf-8"
        )
    )
    gexec.wait()


def generate_sim_params_yaml(
    sim_params_file: str,
    parameters: Dict[str, Any],
    logger: Logger = qclogging.get_basic_logger(),
):
    makedirs(path.dirname(sim_params_file), exist_ok=True)

    logger.debug(f"Raw sim params: {parameters}")
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

    logger.debug(f"Processed sim params: {sim_params}")
    logger.debug(f"Saving sim params to {sim_params_file}")
    with open(sim_params_file, "w") as spf:
        yaml.dump(sim_params, spf)
    logger.debug(f"Sim params saved")


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("realisation_file", type=path.abspath)
    return parser.parse_args()


def main():
    primary_logger = qclogging.get_logger("realisation_to_srf")
    args = load_args()
    rel_df: pd.DataFrame = pd.read_csv(args.realisation_file)
    realisation = rel_df.to_dict(orient="records")[0]
    rel_logger = qclogging.get_realisation_logger(primary_logger, realisation["name"])
    if realisation["type"] == 1:
        create_ps_srf(args.realisation_file, realisation, logger=rel_logger)
    elif realisation["type"] == 2:
        create_ps_ff_srf(args.realisation_file, realisation, logger=rel_logger)
    else:
        raise ValueError(
            f"Type {realisation['type']} faults are not currently supported. "
            f"Contact the software team if you believe this is an error."
        )
    sim_params_file = args.realisation_file.replace(".csv", ".yaml")
    generate_sim_params_yaml(sim_params_file, realisation)


if __name__ == "__main__":
    main()

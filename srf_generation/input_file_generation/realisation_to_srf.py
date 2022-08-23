import argparse
from logging import Logger
from subprocess import run, PIPE, Popen
from typing import Dict, Any, Union, List
from tempfile import NamedTemporaryFile

import yaml
from h5py import File as h5open
import numpy as np
from os import makedirs, path, remove
from qcore import binary_version, srf, geo, qclogging, utils
from qcore.utils import compare_versions
from qcore.uncertainties.mag_scaling import (
    mag2mom,
    MagnitudeScalingRelations,
    mw_to_a_skarlatoudis,
)

from srf_generation.pre_processing_common import (
    calculate_corners,
    get_hypocentre,
    load_realisation_file_as_dict,
)
from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)
from srf_generation.source_parameter_generation.uncertainties.common import (
    HF_RUN_PARAMS,
    BB_RUN_PARAMS,
    LF_RUN_PARAMS,
)


SRF_SUBFAULT_SIZE_KM = 0.1

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


def get_n(fault_size, sub_fault_size):
    return int(round(fault_size / sub_fault_size))


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
    srf2stoch = binary_version.get_unversioned_bin(SRF2STOCH)
    if single_segment:
        command = [srf2stoch, f"target_dx={dx}", f"target_dy={dy}"]
    else:
        command = [srf2stoch, f"dx={dx}", f"dy={dy}"]
    command.extend([f"outfile={stoch_file}", f"infile={srf_file}"])
    logger.debug(f"Creating stoch with command: {command}")
    proc = run(command, stderr=PIPE)
    logger.debug(f"{srf2stoch} stderr: {proc.stderr}")


def get_corners_dbottom(planes, dip_dir=None):
    """
    planes: a list of dictionaries where each dictionary is structured like below.
     {
        "centre": [float(elon), float(elat)],
        "nstrike": int(nstk),
        "ndip": int(ndip),
        "length": float(ln),
        "width": float(wid),
        "strike": stk,
        "dip": dip,
        "shyp": shyp,
        "dhyp": dhyp,
        "dtop": dtop,
    }
    """
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

    return (corners, dbottom)


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

    corners, dbottom = get_corners_dbottom(planes, dip_dir=dip_dir)

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
    tect_type = parameter_dictionary.pop("tect_type", None)
    stype = parameter_dictionary.pop("stype", "cos")
    risetime = parameter_dictionary.pop("risetime", 0.5)
    inittime = parameter_dictionary.pop("inittime", 0.0)

    # srfgen seed is not used by slip2srf (unless we pass rt_rand in),
    # but we also don't need it to be in the sim_params.yaml file
    parameter_dictionary.pop("srfgen_seed", None)

    name = parameter_dictionary.get("name")
    logger.info(f"Generating srf for realisation {name}")

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
        if tect_type is not None and tect_type == "SUBDUCTION_INTERFACE":
            aa = mw_to_a_skarlatoudis(magnitude)
        else:
            # Shallow crustal and subduction slab (currently) use Leonard moment-area scaling relation
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
    with NamedTemporaryFile(mode="w+", delete=False) as gsfp:
        gsfp.write("# nstk= 1 ndip= 1\n")
        gsfp.write(f"# flen= {dd:10.4f} fwid= {dd:10.4f}\n")
        gsfp.write(
            "# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n"
        )
        gsfp.write("1\n")
        gsfp.write(
            f"{longitude:11.5f} {latitude:11.5f} {depth:8.4f} {dd:8.4f} {dd:8.4f} "
            f"{strike:6.1f} {dip:6.1f} {rake:6.1f} {slip:8.2f} {inittime:8.3f}    0\n"
        )

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
    run(commands, stderr=PIPE)

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
        tect_type=tect_type,
    )
    logger.info(
        f"Generated srf for realisation {name}. Moving to next available realisation."
    )


def create_ps_ff_srf(
    realisation_file: str,
    parameter_dictionary: Dict[str, Any],
    stoch_file: Union[None, str] = None,
    logger: Logger = qclogging.get_basic_logger,
):
    name = parameter_dictionary.get("name")
    logger.info(f"Generating srf for realisation {name}")

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
    genslip_version = str(parameter_dictionary.pop("genslip_version"))
    seed = parameter_dictionary.pop("srfgen_seed")
    slip_cov = parameter_dictionary.pop("slip_cov", None)
    risetime_coef = parameter_dictionary.pop("risetime_coef", None)
    rough = parameter_dictionary.pop("rough", 0.0)
    tect_type = parameter_dictionary.pop("tect_type", None)
    vel_mod_1d = parameter_dictionary.pop(
        "srf_vel_mod_1d", DEFAULT_1D_VELOCITY_MODEL_PATH
    )

    asperity_file = parameter_dictionary.pop("asperity_file", None)

    mwsr = MagnitudeScalingRelations(parameter_dictionary.pop("mwsr"))

    # gets
    rvfac = parameter_dictionary.get("rvfac", None)

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
    gen_srf(
        srf_file,
        gsfp.name,
        2,
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
        risetime_coef=risetime_coef,
        rough=rough,
        logger=logger,
        tect_type=tect_type,
        asperity_file=asperity_file,
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
    logger.info(
        f"Generated srf for realisation {name}. Moving to next available realisation."
    )


def create_multi_plane_srf(
    realisation_file: str,
    parameter_dictionary: Dict[str, Any],
    stoch_file: Union[None, str] = None,
    logger: Logger = qclogging.get_basic_logger(),
):
    name = parameter_dictionary.get("name")
    rel_logger = qclogging.get_realisation_logger(logger, name)
    rel_logger.info(f"Generating srf for realisation {name}")

    # pops
    magnitude = parameter_dictionary.pop("magnitude")
    moment = parameter_dictionary.pop("moment")
    rake = parameter_dictionary.pop("rake")
    dip = parameter_dictionary.pop("dip")
    dt = parameter_dictionary.pop("dt", 0.005)
    dtop = parameter_dictionary.pop("dtop")
    length = parameter_dictionary.pop("length")
    shypo = parameter_dictionary.pop("shypo")
    dhypo = parameter_dictionary.pop("dhypo")
    genslip_version = str(parameter_dictionary.pop("genslip_version"))
    dip_dir = parameter_dictionary.pop("dip_dir")
    seed = parameter_dictionary.pop("srfgen_seed")
    tect_type = parameter_dictionary.pop("tect_type")
    fault_type = parameter_dictionary.pop("fault_type")
    plane_count = parameter_dictionary.pop("plane_count")

    rough = parameter_dictionary.pop("rough", 0.0)
    slip_cov = parameter_dictionary.pop("slip_cov", None)

    asperity_file = parameter_dictionary.pop("asperity_file", None)

    strike = [
        parameter_dictionary.pop(f"strike_subfault_{i}") for i in range(plane_count)
    ]
    flen = [
        parameter_dictionary.pop(f"length_subfault_{i}") for i in range(plane_count)
    ]
    fwid = [parameter_dictionary.pop(f"width_subfault_{i}") for i in range(plane_count)]

    clon = [parameter_dictionary.pop(f"clon_subfault_{i}") for i in range(plane_count)]
    clat = [parameter_dictionary.pop(f"clat_subfault_{i}") for i in range(plane_count)]

    for key in list(parameter_dictionary.keys()):
        # Remove all the subfault keys so they don't get passed to the sim_params.yaml file
        # Must make a copy of the keys to be able to remove them.
        if "_subfault_" in key:
            parameter_dictionary.pop(key)

    # gets
    vel_mod_1d = parameter_dictionary.pop(
        "srf_vel_mod_1d", DEFAULT_1D_VELOCITY_MODEL_PATH
    )
    rvfac = parameter_dictionary.get("rvfac", None)

    logger.debug(
        "All srf generation parameters successfully obtained from the realisation file"
    )

    dlen = dwid = SRF_SUBFAULT_SIZE_KM

    nx = [get_n(flen[i], dlen) for i in range(plane_count)]
    ny = get_n(fwid[0], dwid)

    if (
        utils.compare_versions(genslip_version, "5.4.2") < 0
        and tect_type == "SUBDUCTION_INTERFACE"
    ):
        raise RuntimeError(
            "Subduction interface faults are only available for version 5.4.2 and above"
        )

    gsf_file = NamedTemporaryFile(mode="w", delete=False)
    rel_logger.debug(f"Gsf will be saved to the temporary file {gsf_file.name}")
    srf_file = realisation_file.replace(".csv", ".srf")
    rel_logger.debug(f"Srf will be saved to {srf_file}")
    corners_file = realisation_file.replace(".csv", ".corners")
    rel_logger.debug(f"The corners file will be saved to {corners_file}")

    with NamedTemporaryFile(mode="w", delete=False) as gsfp:
        rel_logger.debug("Saving segments file to {}".format(gsfp.name))
        gsfp.write("{}\n".format(plane_count))
        for f in range(plane_count):
            gsfp.write(
                "{:f} {:f} {:f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:d} {:d}\n".format(
                    clon[f],
                    clat[f],
                    dtop,
                    strike[f],
                    dip,
                    rake,
                    flen[f],
                    fwid[f],
                    nx[f],
                    ny,
                )
            )

    fault_seg_bin = binary_version.get_unversioned_bin(FAULTSEG2GSFDIPDIR)
    cmd = [
        fault_seg_bin,
        "read_slip_vals=0",
        "infile={}".format(gsfp.name),
        "outfile={}".format(gsf_file.name),
    ]
    if dip_dir is not None:
        cmd.append("dipdir={}".format(dip_dir))

    rel_logger.info(
        "Calling fault_seg2gsf_dipdir with command {}".format(" ".join(cmd))
    )
    gexec = run(cmd)
    rel_logger.debug(f"{fault_seg_bin} finished running with stderr: {gexec.stderr}")

    # remove(gsfp.name)
    # rel_logger.debug("Removed segments file")

    if int(plane_count > 1):
        rel_logger.debug("Multiple segments detected. Generating xseg argument")
        flen_array = np.asarray(flen)
        xseg = ",".join(map(str, flen_array.cumsum() - flen_array / 2))
    else:
        xseg = "-1"

    gen_srf(
        srf_file,
        gsf_file.name,
        4,
        magnitude,
        dt,
        sum(nx),
        ny,
        seed,
        shypo,
        dhypo,
        vel_mod_1d,
        genslip_version=genslip_version,
        rvfac=rvfac,
        rough=rough,
        slip_cov=slip_cov,
        xseg=xseg,
        logger=rel_logger,
        tect_type=tect_type,
        asperity_file=asperity_file,
    )

    rel_logger.info("srf generated, creating stoch")

    if stoch_file is None:
        stoch_file = realisation_file.replace(".csv", ".stoch")
    create_stoch(stoch_file, srf_file, single_segment=(plane_count == 1), logger=logger)

    rel_logger.info("stoch created, making info")

    # save INFO
    create_info_file(
        srf_file,
        4,
        magnitude,
        rake,
        dt,
        tect_type=tect_type,
        dip_dir=dip_dir,
        shypo=[shypo],
        dhypo=dhypo,
        vm=vel_mod_1d,
        logger=rel_logger,
    )

    logger.info(
        f"Generated srf for realisation {name}. Moving to next available realisation."
    )

    # path to resulting SRF
    return srf_file


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
    type,
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
    rough=0.0,
    slip_cov=None,
    risetime_coef=None,
    tect_type=None,
    fault_planes=1,
    asperity_file=None,
    xseg: Union[float, List[float]] = "-1",
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    :param xseg: Genslip parameter:
        For multi plane arrays the length (along strike) of each plane. -1 for single plane.
        Deprecated in genslip 5.4.2, ignored if present
    """
    genslip_bin = binary_version.get_genslip_bin(genslip_version)
    if (
        compare_versions(genslip_version, "5.4.2") < 0
        and tect_type == "SUBDUCTION_INTERFACE"
    ):
        raise AssertionError(
            "Cannot generate subduction srfs with genslip version less than 5.4.2"
        )
    if compare_versions(genslip_version, "5") > 0:
        # Positive so version greater than 5
        logger.debug(
            "Using genslip version {}. Using nstk and ndip and rup_delay (for type 4)".format(
                genslip_version
            )
        )
        xstk = "nstk"
        ydip = "ndip"
        rup_name = "rup_delay"
    else:
        # Not positive so version at most 5
        logger.debug(
            "Using genslip version {}. Using nx and ny and rupture_delay (for type 4)".format(
                genslip_version
            )
        )
        xstk = "nx"
        ydip = "ny"
        rup_name = "rupture_delay"
    cmd = [
        genslip_bin,
        "read_erf=0",
        "write_srf=1",
        "read_gsf=1",
        "write_gsf=0",
        f"infile={gsf_file}",
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
    if type == 4:
        cmd.extend(
            [
                f"seg_delay={0}",
                f"nseg={fault_planes}",
                f"nseg_bounds={fault_planes - 1}",
                f"xseg={xseg}",
                f"rvfac_seg=-1",
                f"gwid=-1",
                f"side_taper=0.02",
                f"bot_taper=0.02",
                f"top_taper=0.0",
                f"{rup_name}=0",
            ]
        )

    if tect_type == "SUBDUCTION_INTERFACE":
        cmd.extend(
            [
                "kmodel=-1",
                "xmag_exp=0.5",
                "ymag_exp=0.5",
                "kx_corner=2.5482",
                "ky_corner=2.3882",
                "tsfac_slope=-0.5",
                "tsfac_bzero=-0.1",
            ]
        )
        if risetime_coef is None:
            cmd.append(f"risetime_coef=1.95")
    if rvfac is not None:
        cmd.append(f"rvfrac={rvfac}")
    if rough is not None:
        cmd.append(f"alpha_rough={rough}")
    if slip_cov is not None:
        cmd.append(f"slip_sigma={slip_cov}")
    if risetime_coef is not None:
        cmd.append(f"risetime_coef={risetime_coef}")
    if asperity_file is not None:
        cmd.append(f"read_slip_file=1")
        cmd.append(f"init_slip_file={asperity_file}")
    logger.debug("Creating SRF with command: {}".format(" ".join(cmd)))
    with open(srf_file, "w") as srfp:
        proc = run(cmd, stdout=srfp, stderr=PIPE)
    logger.debug(f"{genslip_bin} stderr: {proc.stderr}")


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
        f"1\n{lon:f} {lat:f} {dtop:f} {strike} {dip} {rake} {flen:f} {fwid:f} {nx:d} {ny:d}".encode(
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
        elif key == "vs30_file_path":
            sim_params["stat_vs_est"] = value
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
    qclogging.add_general_file_handler(primary_logger, "rel2srf.txt")
    args = load_args()
    realisation = load_realisation_file_as_dict(args.realisation_file)
    rel_logger = qclogging.get_realisation_logger(primary_logger, realisation["name"])
    if realisation["type"] == 1:
        create_ps_srf(args.realisation_file, realisation, logger=rel_logger)
    elif realisation["type"] == 2:
        create_ps_ff_srf(args.realisation_file, realisation, logger=rel_logger)
    elif realisation["type"] == 4:
        create_multi_plane_srf(args.realisation_file, realisation, logger=rel_logger)
    else:
        raise ValueError(
            f"Type {realisation['type']} faults are not currently supported. "
            f"Contact the software team if you believe this is an error."
        )
    sim_params_file = args.realisation_file.replace(".csv", ".yaml")
    generate_sim_params_yaml(sim_params_file, realisation)


if __name__ == "__main__":
    main()

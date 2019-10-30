import argparse
from subprocess import call, PIPE
from h5py import File as h5open
import numpy as np
from os import makedirs, path
import pandas as pd
from qcore import binary_version, srf, geo, simulation_structure

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


def gen_meta(
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
        file_name = path.splitext(srf_file)[0]
    file_name = "{}.info".format(file_name)
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
    name,
    latitude,
    longitude,
    depth,
    magnitude,
    strike,
    rake,
    dip,
    cs_root,
    moment=0,
    dt=0.005,
    stoch=None,
    vs=3.20,
    rho=2.44,
    target_area_km=None,
    target_slip_cm=None,
    stype="cos",
    rise_time=0.5,
    init_time=0.0,
    silent=False,
    **kwargs,
):
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
        aa = np.exp(2.0 * np.log(moment) / 3.0 - 14.7 * np.log(10.0))
        dd = np.sqrt(aa)
        slip = (moment * 1.0e-20) / (aa * vs * vs * rho)

    ###
    ### file names
    ###
    srf_file = simulation_structure.get_srf_path(cs_root, name)
    gsf_file = srf_file.replace(".srf", ".gsf")
    makedirs(path.dirname(srf_file), exist_ok=True)

    ###
    ### create GSF
    ###
    with open(gsf_file, "w") as gsfp:
        gsfp.write("# nstk= 1 ndip= 1\n")
        gsfp.write(f"# flen= {dd:10.4f} fwid= {dd:10.4f}\n")
        gsfp.write(
            "# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n"
        )
        gsfp.write("1\n")
        gsfp.write(
            f"{longitude:11.5f} {latitude:11.5f} {depth:8.4f} {dd:8.4f} {dd:8.4f} "
            f"{strike:6.1f} {dip:6.1f} {rake:6.1f} {slip:8.2} {init_time:8.3f}    0\n"
        )

    ###
    ### create SRF
    ###
    if silent:
        stderr = PIPE
    else:
        stderr = None
    commands = [
        binary_version.get_unversioned_bin(GENERICSLIP2SRF),
        f"infile={gsf_file}",
        f"outfile={srf_file}",
        "outbin=0",
        f"stype={stype}",
        f"dt={dt}",
        "plane_header=1",
        f"risetime={rise_time}",
        "risetimefac=1.0",
        "risetimedep=0.0",
    ]
    call(commands, stderr=stderr)

    ###
    ### save STOCH
    ###
    stoch_file = simulation_structure.get_stoch_path(cs_root, name)
    gen_stoch(stoch_file, srf_file)

    ###
    ### save INFO
    ###
    gen_meta(
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

    # location of resulting SRF file
    return srf_file


def load_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("srfgenparams_file", type=path.abspath)
    parser.add_argument(
        "-c", "--cybershake_root", type=path.abspath, default=path.abspath(".")
    )
    return parser.parse_args()


def main():
    args = load_args()
    rel_df: pd.DataFrame = pd.read_csv(args.srfgenparams_file)
    realisation = rel_df.to_dict(orient='records')[0]
    if realisation["type"] == 1:
        create_ps_srf(cs_root=args.cybershake_root, **realisation)


if __name__ == "__main__":
    main()

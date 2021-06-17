from argparse import ArgumentParser
from logging import Logger
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import yaml

from qcore import geo, gmt, qclogging


def plot_vm(
    vm_params_dict: dict,
    srf_corners: np.ndarray,
    land_outline_path: Path,
    centre_line_path: Path,
    mag: float,
    outdir: Path,
    ptemp: Path,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Plots VM domain as well as SRF domain if possible

    Parameters
    ----------
    vm_params_dict :
    srf_corners :
    land_outline_path :
    centre_line_path :
    mag :
    outdir :
    ptemp :
    logger :
    """

    logger.debug("Plotting vm")
    p = gmt.GMTPlot(ptemp / "optimisation.ps")
    p.spacial("M", vm_params_dict["plot_region"], sizing=7)
    p.coastlines()

    srf_path = ptemp / "srf.path"
    if srf_path.exists():
        # filled slip area
        p.path(srf_path, is_file=True, fill="yellow", split="-")
        # top edge
        for plane in srf_corners:
            p.path(
                "\n".join([" ".join(map(str, ll)) for ll in plane[:2]]), is_file=False
            )

    # vm domain (simple and adjusted)
    p.path(vm_params_dict["path"], is_file=False, close=True, fill="black@95")
    if vm_params_dict["adjusted"]:
        p.path(
            vm_params_dict["path_mod"],
            is_file=False,
            close=True,
            fill="black@95",
            split="-",
            width="1.0p",
        )

    # info text for simple and adjusted domains
    p.text(
        sum(vm_params_dict["plot_region"][0:2]) / 2.0,
        vm_params_dict["plot_region"][3],
        "Mw: %.2f X: %.0fkm, Y: %.0fkm, land: %.0f%%"
        % (mag, vm_params_dict["xlen"], vm_params_dict["ylen"], vm_params_dict["land"]),
        align="CT",
        dy=-0.1,
        box_fill="white@50",
    )
    if vm_params_dict["adjusted"]:
        p.text(
            sum(vm_params_dict["plot_region"][0:2]) / 2.0,
            vm_params_dict["plot_region"][3],
            "MODIFIED land: %.0f%%" % (vm_params_dict["land_mod"]),
            align="CT",
            dy=-0.25,
            box_fill="white@50",
        )

    # land outlines blue, nz centre line (for bearing calculation) red
    p.path(land_outline_path, is_file=True, close=False, colour="blue", width="0.2p")
    p.path(centre_line_path, is_file=True, close=False, colour="red", width="0.2p")

    # actual corners retrieved from NZVM output or generated if args.novm
    # not available if VM was skipped
    corner_file = outdir / "VeloCorners.txt"

    if corner_file.exists():
        logger.debug("Getting corners from VeloModCorners.txt")
        p.points(
            corner_file,
            fill="red",
            line=None,
            shape="c",
            size=0.05,
        )
    else:
        logger.debug("VeloModCorners.txt doesn't exist, deriving corners from path mod")
        p.points(
            vm_params_dict["path_mod"],
            is_file=False,
            fill="red",
            line=None,
            shape="c",
            size="0.05",
        )

    # store PNG
    p.finalise()
    logger.debug("Saving image")

    p.png(
        dpi=200,
        clip=True,
        background="white",
        out_name=(outdir / vm_params_dict["name"]).resolve(),
    )


def main(
    name: str,
    vm_params_dict: dict,
    outdir: Path,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    vm_params_dict loaded from vm_params.yaml doesn't have all info plot_vm() needs.
    This function gathers and works out the necessary input (except SRF-relevant info) to run this file as a stand-alone script

    Parameters
    ----------
    name : name of the fault/event
    vm_params_dict : Dictionary extracted from vm_params.yaml
    outdir :
    logger :
    """
    from rel2vm_params import (
        get_vm_land_proportion,
        corners2region,
        NZ_CENTRE_LINE,
        NZ_LAND_OUTLINE,
    )

    # vm_params_dict is the dictionary directly loaded from vm_params.yaml

    with TemporaryDirectory(prefix=f"_tmp_{name}_", dir=outdir) as ptemp:
        ptemp = Path(ptemp)
        vm_params_dict["name"] = name

        origin = (vm_params_dict["MODEL_LON"], vm_params_dict["MODEL_LAT"])
        xlen = vm_params_dict["extent_x"]
        ylen = vm_params_dict["extent_y"]
        c1, c2, c3, c4 = geo.build_corners(
            origin, vm_params_dict["MODEL_ROT"], xlen, ylen
        )

        vm_params_dict[
            "path_mod"
        ] = "{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n".format(
            c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]
        )
        vm_params_dict["path"] = vm_params_dict["path_mod"]
        vm_params_dict["adjusted"] = False
        vm_params_dict["xlen"] = xlen
        vm_params_dict["ylen"] = ylen
        vm_params_dict["land"] = get_vm_land_proportion(c1, c2, c3, c4)

        vm0_region = corners2region(c1, c2, c3, c4)
        plot_region = (
            vm0_region[0] - 1,
            vm0_region[1] + 1,
            vm0_region[2] - 1,
            vm0_region[3] + 1,
        )

        vm_params_dict["plot_region"] = plot_region

        # plotting the domain of VM. No SRF
        plot_vm(
            vm_params_dict,
            [],
            NZ_LAND_OUTLINE,
            NZ_CENTRE_LINE,
            vm_params_dict["mag"],
            outdir,
            ptemp,
            logger=logger,
        )


def load_args(logger: Logger = qclogging.get_basic_logger()):
    """
    Unpacks arguments and does basic checks

    Parameters
    ----------
    logger :

    Returns
    -------

    """
    parser = ArgumentParser()
    arg = parser.add_argument

    arg("name", help="Name of the fault")

    arg(
        "vm_params_path",
        help="path to vm_params.yaml",
    )

    arg(
        "-o",
        "--outdir",
        help="output directory"
        "(if not specified, the same location as vm_params.yaml is in",
        default=None,
    )

    args = parser.parse_args()
    args.vm_params_path = Path(args.vm_params_path).resolve()

    if args.outdir is None:
        args.outdir = args.vm_params_path.parent
    args.outdir = Path(args.outdir).resolve()

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("plot_vm")
    qclogging.add_general_file_handler(logger, Path.cwd() / "plot_vm.txt")
    args = load_args(logger=logger)

    with open(args.vm_params_path, "r") as f:
        vm_params_dict = yaml.load(f, Loader=yaml.SafeLoader)

    main(args.name, vm_params_dict, args.outdir, logger=logger)

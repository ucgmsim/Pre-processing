"""
Plots the VM domain and the SRF planes
"""

from argparse import ArgumentParser
from io import StringIO
from logging import Logger
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import yaml

from qcore import geo, gmt, qclogging, validate_vm


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
    vm_params_dict : dict
        Dictionary extracted from vm_params.yaml file.
    srf_corners : np.ndarray.
        Corners of SRF planes. Formatted as [[[lon,lat],[lon,lat],[lon,lat],[lon,lat]],...]
    land_outline_path : Path
        Path to the land outline file
    centre_line_path : Path
        Path to the centre line file
    mag : float
        Magnitude of the event
    outdir : Path
        Output directory
    ptemp : Path
        Temporary directory
    logger : Logger
        Default is qclogging.get_basic_logger()
    """

    from rel2vm_params import write_srf_path

    logger.debug("Plotting vm")
    p = gmt.GMTPlot(ptemp / "optimisation.ps")
    p.spacial("M", vm_params_dict["plot_region"], sizing=7)
    p.coastlines()

    # SRF domain
    if len(srf_corners) > 0:
        srf_path = write_srf_path(
            srf_corners, ptemp
        )  # write srf.path file if not already present, otherwise use the existing file

        # filled slip area
        p.path(srf_path, is_file=True, fill="yellow", split="-")
        # top edge
        for plane in srf_corners:
            p.path(
                "\n".join([" ".join(map(str, ll)) for ll in plane[:2]]), is_file=False
            )

    # plot the VM domain (simple and adjusted)
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
        p.points(corner_file, fill="red", line=None, shape="c", size=0.05)
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
    rel_path: Path,
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Gathers necessary input to call plot_vm() function to plot VM domain and SRF planes (if realisation CSV is supplied)

    Parameters
    ----------
    name : str
        name of the fault/event. This is used to name the output file.
    vm_params_dict : dict
        Dictionary extracted from vm_params.yaml file.
    outdir : Path
        Output directory
    rel_path : Path
        Path to the realisation csv file.

    logger :
    """
    from rel2vm_params import (
        get_vm_land_proportion,
        corners2region,
        NZ_CENTRE_LINE,
        NZ_LAND_OUTLINE,
        load_rel,
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

        vm_params_dict["path_mod"] = (
            "{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n{:.6f}\t{:.6f}\n".format(
                c1[0], c1[1], c2[0], c2[1], c3[0], c3[1], c4[0], c4[1]
            )
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

        srf_meta = load_rel(rel_path)
        srf_corners = srf_meta["corners"]

        # plotting the domain of VM.
        plot_vm(
            vm_params_dict,
            srf_corners,
            NZ_LAND_OUTLINE,
            NZ_CENTRE_LINE,
            vm_params_dict["mag"],
            outdir,
            ptemp,
            logger=logger,
        )

        # Validate the VM domain. Only needed when this code is run as a standalone script.
        # If this code is run as part of the VM workflow, the validation is done in the main script.
        polygon = np.loadtxt(StringIO(vm_params_dict["path"]))
        errors = validate_vm.validate_region(polygon)
        errors.extend(validate_vm.validate_vm_bounds(polygon, srf_corners))
        if errors:
            logger.warning(f"WARNING: {errors}")


def load_args(logger: Logger = qclogging.get_basic_logger()):
    """
    Unpacks arguments and does basic checks

    Parameters
    ----------
    logger :  Logger
        Default is qclogging.get_basic_logger()

    Returns
    -------
    Processed arguments

    """
    parser = ArgumentParser()
    arg = parser.add_argument

    arg("vm_params_path", help="path to vm_params.yaml", type=Path)
    arg("rel_path", help="Path to the realisation csv file", type=Path)

    arg(
        "-o",
        "--outdir",
        help="output directory"
        "(if not specified, the same location as vm_params.yaml is in",
        default=None,
    )

    args = parser.parse_args()
    args.vm_params_path = args.vm_params_path.resolve()
    args.rel_path = args.rel_path.resolve()

    if (
        args.outdir is None
    ):  # if not specified, use the directory that contains vm_params.yaml
        args.outdir = args.vm_params_path.parent

    args.outdir = Path(args.outdir).resolve()

    args.outdir.mkdir(exist_ok=True, parents=True)
    assert args.vm_params_path.exists(), f"File is not present: {args.vm_params_path}"
    assert args.rel_path.exists(), f"File is not present: {args.rel_path}"

    args.name = args.rel_path.stem

    return args


if __name__ == "__main__":
    logger = qclogging.get_logger("plot_vm")
    qclogging.add_general_file_handler(logger, Path.cwd() / "plot_vm.txt")
    args = load_args(logger=logger)

    with open(args.vm_params_path, "r") as f:
        vm_params_dict = yaml.load(f, Loader=yaml.SafeLoader)

    main(args.name, vm_params_dict, args.outdir, args.rel_path, logger=logger)

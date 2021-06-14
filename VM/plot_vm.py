import os
from qcore import gmt, qclogging
from logging import Logger


def plot_vm(
    vm_params,
    srf_corners,
    land_outline,
    centre_line,
    mag,
    outdir,
    ptemp,
    logger: Logger = qclogging.get_basic_logger(),
):
    """

    Parameters
    ----------
    vm_params :
    srf_corners :
    mag :
    outdir :
    ptemp :
    logger :
    """

    logger.debug("Plotting vm")
    p = gmt.GMTPlot(os.path.join(ptemp, "optimisation.ps"))
    p.spacial("M", vm_params["plot_region"], sizing=7)
    p.coastlines()

    # filled slip area
    p.path("%s/srf.path" % (ptemp), is_file=True, fill="yellow", split="-")
    # top edge
    for plane in srf_corners:
        p.path("\n".join([" ".join(map(str, ll)) for ll in plane[:2]]), is_file=False)

    # vm domain (simple and adjusted)
    p.path(vm_params["path"], is_file=False, close=True, fill="black@95")
    if vm_params["adjusted"]:
        p.path(
            vm_params["path_mod"],
            is_file=False,
            close=True,
            fill="black@95",
            split="-",
            width="1.0p",
        )

    # info text for simple and adjusted domains
    p.text(
        sum(vm_params["plot_region"][0:2]) / 2.0,
        vm_params["plot_region"][3],
        "Mw: %.2f X: %.0fkm, Y: %.0fkm, land: %.0f%%"
        % (mag, vm_params["xlen"], vm_params["ylen"], vm_params["land"]),
        align="CT",
        dy=-0.1,
        box_fill="white@50",
    )
    if vm_params["adjusted"]:
        p.text(
            sum(vm_params["plot_region"][0:2]) / 2.0,
            vm_params["plot_region"][3],
            "MODIFIED land: %.0f%%" % (vm_params["land_mod"]),
            align="CT",
            dy=-0.25,
            box_fill="white@50",
        )

    # land outlines blue, nz centre line (for bearing calculation) red
    p.path(land_outline, is_file=True, close=False, colour="blue", width="0.2p")
    p.path(centre_line, is_file=True, close=False, colour="red", width="0.2p")

    # actual corners retrieved from NZVM output or generated if args.novm
    # not available if VM was skipped
    vm_exists = "vm_dir" in vm_params
    if vm_exists:
        logger.debug("vm exists, getting corners from VeloModCorners")
        p.points(
            "%s/VeloModCorners.txt" % (vm_params["vm_dir"]),
            fill="red",
            line=None,
            shape="c",
            size=0.05,
        )
    else:
        logger.debug("vm doesn't exist, deriving corners from path mod")
        p.points(
            vm_params["path_mod"],
            is_file=False,
            fill="red",
            line=None,
            shape="c",
            size="0.05",
        )

    # store PNG
    p.finalise()
    logger.debug("Saving image")
    if vm_exists:
        p.png(dpi=200, clip=True, background="white", out_dir=vm_params["vm_dir"])
    else:
        p.png(
            dpi=200,
            clip=True,
            background="white",
            out_name=os.path.abspath(os.path.join(outdir, vm_params["name"])),
        )

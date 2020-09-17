#!/usr/bin/env python2

import os
from logging import Logger
from subprocess import check_call
import sys

from qcore import qclogging
from qcore.binary_version import get_unversioned_bin
from qcore.constants import VM_PARAMS_FILE_NAME
from qcore.utils import load_yaml


def gen_coords(
    vm_dir=".",
    debug=False,
    geoproj="1",
    do_coords="1",
    centre_origin="1",
    logger: Logger = qclogging.get_basic_logger(),
):
    """
    Generate coordinate files for an emod3d domain (set of grid points).
    outdir: directory to store coordinate files in
    debug: print additional info
    geoproj:
    do_coords:
    """

    logger.debug("Generating coords. VM directory: {}".format(vm_dir))

    # load params for velocity model
    vm = load_yaml(os.path.join(vm_dir, VM_PARAMS_FILE_NAME))

    XLEN = vm["nx"] * vm["hh"]
    YLEN = vm["ny"] * vm["hh"]
    ZLEN = vm["nz"] * vm["hh"]

    # list of outputs that this function can create
    GRIDFILE = os.path.join(vm_dir, "gridfile%s" % (vm["sufx"]))
    GRIDOUT = os.path.join(vm_dir, "gridout%s" % (vm["sufx"]))
    MODEL_COORDS = os.path.join(vm_dir, "model_coords%s" % (vm["sufx"]))
    MODEL_PARAMS = os.path.join(vm_dir, "model_params%s" % (vm["sufx"]))
    MODEL_BOUNDS = os.path.join(vm_dir, "model_bounds%s" % (vm["sufx"]))

    # generate gridfile
    try:
        with open(GRIDFILE, "w") as gridf:
            gridf.write("xlen=%f\n" % (XLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, XLEN, vm["hh"]))
            gridf.write("ylen=%f\n" % (YLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, YLEN, vm["hh"]))
            gridf.write("zlen=%f\n" % (ZLEN))
            gridf.write("%10.4f %10.4f %13.6e\n" % (0.0, ZLEN, vm["hh"]))
    except IOError:
        message = "Cannot write GRIDFILE: {}".format(GRIDFILE)
        logger.log(qclogging.NOPRINTERROR, message)
        raise IOError(message)
    logger.debug("Wrote grid file to {}".format(GRIDFILE))

    # generate model_params
    cmd = (
        "{} "
        "geoproj={geoproj} gridfile='{GRIDFILE}' gridout='{GRIDOUT}' "
        "center_origin={centre_origin} do_coords={do_coords} "
        "nzout=1 name='{MODEL_COORDS}' gzip=0 latfirst=0 "
        "modellon={vm[MODEL_LON]} modellat={vm[MODEL_LAT]} "
        "modelrot={vm[MODEL_ROT]} 1> '{MODEL_PARAMS}'"
    ).format(get_unversioned_bin("gen_model_cords"), **dict(locals(), **globals()))
    logger.log(qclogging.logging.DEBUG * (1 + int(debug)), cmd)
    if not debug:
        cmd += " 2>/dev/null"
    check_call(cmd, shell=True)

    # also generate coordinate related outputs
    if do_coords != "1":
        return

    # retrieve MODEL_BOUNDS
    x_bounds = [0, vm["nx"] - 1]
    y_bounds = [0, vm["ny"] - 1]
    try:
        with open(MODEL_COORDS, "r") as coordf:
            with open(MODEL_BOUNDS, "w") as boundf:
                for line in coordf:
                    x, y = map(float, line.split()[2:4])
                    if x in x_bounds or y in y_bounds:
                        boundf.write(line)
    except IOError:
        message = "Cannot write MODEL_BOUNDS: {}".format(MODEL_BOUNDS)
        logger.log(qclogging.NOPRINTERROR, message)
        raise IOError(message)


# allow running from shell
if __name__ == "__main__":
    if len(sys.argv) > 1:
        gen_coords(vm_dir=sys.argv[1])
    else:
        gen_coords()

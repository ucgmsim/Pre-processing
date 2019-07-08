#!/usr/bin/env python3

import os
from subprocess import check_call
import sys

from qcore.binary_version import get_unversioned_bin
from qcore.constants import VM_PARAMS_FILE_NAME
from qcore.utils import load_yaml


def gen_coords(vm_dir=".", debug=False, geoproj="1", do_coords="1", centre_origin="1"):
    """
    Generate coordinate files for an emod3d domain (set of grid points).
    outdir: directory to store coordinate files in
    debug: print additional info
    geoproj: 
    do_coords: 
    """

    # load params for velocity model
    vm = load_yaml(os.path.join(vm_dir, VM_PARAMS_FILE_NAME))

    hh = vm["hh"]
    XLEN = vm["nx"] * hh
    YLEN = vm["ny"] * hh
    ZLEN = vm["nz"] * hh

    # list of outputs that this function can create
    GRIDFILE = os.path.join(vm_dir, "gridfile{}".format(vm["sufx"]))
    GRIDOUT = os.path.join(vm_dir, "gridout{}".format(vm["sufx"]))
    MODEL_COORDS = os.path.join(vm_dir, "model_coords{}".format(vm["sufx"]))
    MODEL_PARAMS = os.path.join(vm_dir, "model_params{}".format(vm["sufx"]))
    MODEL_BOUNDS = os.path.join(vm_dir, "model_bounds{}".format(vm["sufx"]))

    # generate gridfile
    try:
        with open(GRIDFILE, "w") as gridf:
            gridf.write("xlen={:f}\n".format(XLEN))
            gridf.write("{:10.4f} {:10.4f} {:13.6e}\n".format(0.0, XLEN, hh))
            gridf.write("ylen={:f}\n".format(YLEN))
            gridf.write("{:10.4f} {:10.4f} {:13.6e}\n".format(0.0, YLEN, hh))
            gridf.write("zlen={:f}\n".format(ZLEN))
            gridf.write("{:10.4f} {:10.4f} {:13.6e}\n".format(0.0, ZLEN, hh))
    except IOError:
        raise IOError("Cannot write GRIDFILE: %s" % (GRIDFILE))

    # generate model_params
    cmd = (
        "{} "
        "geoproj={geoproj} gridfile='{GRIDFILE}' gridout='{GRIDOUT}' "
        "center_origin={centre_origin} do_coords={do_coords} "
        "nzout=1 name='{MODEL_COORDS}' gzip=0 latfirst=0 "
        "modellon={vm[MODEL_LON]} modellat={vm[MODEL_LAT]} "
        "modelrot={vm[MODEL_ROT]} 1> '{MODEL_PARAMS}'"
    ).format(get_unversioned_bin("gen_model_cords"), **locals())
    if debug:
        print(cmd)
    else:
        cmd += " 2>/dev/null"
    check_call(cmd, shell=True)

    # also generate coordinate related outputs
    if do_coords == "1":

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
            raise IOError("Cannot write MODEL_BOUNDS: {}".format(MODEL_BOUNDS))


# allow running from shell
if __name__ == "__main__":
    if len(sys.argv) > 1:
        gen_coords(vm_dir=sys.argv[1])
    else:
        gen_coords()

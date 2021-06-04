import csv
from logging import Logger
import os

from qcore import qclogging

script_dir = os.path.dirname(os.path.abspath(__file__))
NZ_CENTRE_LINE = os.path.join(script_dir, "../SrfGen/NHM/res/centre.txt")
NZ_LAND_OUTLINE = os.path.join(script_dir, "../SrfGen/NHM/res/rough_land.txt")


def temp_paths(ptemp):
    vm_working_dir = os.path.join(ptemp, "output")
    nzvm_cfg = os.path.join(ptemp, "nzvm.cfg")
    vm_params_path = os.path.join(ptemp, "vm_params.yaml")
    return (vm_working_dir, nzvm_cfg, vm_params_path)


def store_summary(table, info_store, logger: Logger = qclogging.get_basic_logger()):
    # initialise table file
    logger.debug("Saving summary")
    keys = info_store[0].keys()
    with open(table, "w") as csvfile:
        writer = csv.writer(csvfile, delimiter=",")
        writer.writerow(keys)
        for i in info_store:
            writer.writerow([i[key] for key in keys])
    logger.info("Saved summary to {}".format(table))

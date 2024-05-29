#!/usr/bin/env python3
import collections
import dataclasses
import re
from pathlib import Path
from typing import Generator, Tuple

import numpy as np
import qcore.coordinates
import scipy as sp
import yaml
from nshmdb import fault
from nshmdb.fault import Fault

from srf_generation.source_parameter_generation.common import \
    DEFAULT_1D_VELOCITY_MODEL_PATH


@dataclasses.dataclass
class FaultJump:
    jump_location_lat: float
    jump_location_lon: float


@dataclasses.dataclass
class RealisationFault(Fault):
    magnitude: float = None
    parent: "RealisationFault" = None
    parent_jump_coords: Tuple[float, float] | None = None
    shyp: float = None
    dhyp: float = None

    def random_fault_coordinates(self) -> (float, float):
        weibull_scale = 0.612
        dhyp = (
            self.widths()[0]
            * sp.stats.weibull_min(3.353, 0, 1, scale=weibull_scale).rvs(1)[0]
        )
        shyp = distributions.rand_shyp() * np.sum(self.lengths())
        return (shyp, dhyp)

    def expected_fault_coordinates(self) -> (float, float):
        weibull_scale = 0.612
        dhyp = (
            self.widths()[0]
            * sp.stats.truncweibull_min(3.353, 0, 1, scale=weibull_scale).expect()
        )
        shyp = 0
        return (shyp, float(dhyp))


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    dt: float
    magnitude: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    velocity_model: str
    faults: dict[str, RealisationFault]

    def initial_fault(self) -> RealisationFault:
        return next(fault for fault in self.faults.values() if not fault.parent)


def read_realisation(realisation_filepath: Path) -> Realisation:
    with open(realisation_filepath, "r", encoding="utf-8") as realisation_file:
        raw_yaml_data = yaml.safe_load(realisation_file)
        faults_object = raw_yaml_data["faults"]
        faults = {
            name: RealisationFault(
                name=fault_obj["name"],
                tect_type=fault_obj["tect_type"],
                shyp=fault_obj["shyp"],
                dhyp=fault_obj["dhyp"],
                magnitude=fault_obj["magnitude"],
                parent_jump_coords=(
                    tuple(fault_obj["parent_jump_coords"])
                    if fault_obj["parent_jump_coords"]
                    else None
                ),
                planes=[
                    fault.FaultPlane(
                        qcore.coordinates.wgs_depth_to_nztm(
                            np.array(params["corners"])
                        ),
                        params["rake"],
                    )
                    for params in fault_obj["planes"]
                ],
            )
            for name, fault_obj in faults_object.items()
        }

        for name, fault_obj in faults_object.items():
            if fault_obj["parent"]:
                faults[name].parent = faults[fault_obj["parent"]]

        return Realisation(
            name=raw_yaml_data["name"],
            type=raw_yaml_data["type"],
            magnitude=fault_obj["magnitude"],
            dt=raw_yaml_data["dt"],
            genslip_seed=raw_yaml_data["genslip_seed"],
            srfgen_seed=raw_yaml_data["srfgen_seed"],
            genslip_version=raw_yaml_data["genslip_version"],
            faults=faults,
            velocity_model=DEFAULT_1D_VELOCITY_MODEL_PATH,
        )


def topologically_sorted_faults(
    faults: list[RealisationFault],
) -> Generator[RealisationFault, None, None]:
    fault_children_map = collections.defaultdict(list)
    for fault in faults:
        if fault.parent:
            fault_children_map[fault.parent.name].append(fault)

    def bfs_traverse_fault_map(
        fault: RealisationFault,
    ) -> Generator[RealisationFault, None, None]:
        yield fault
        for child in fault_children_map[fault.name]:
            yield from bfs_traverse_fault_map(child)

    initial_fault = next(fault for fault in faults if not fault.parent)
    yield from bfs_traverse_fault_map(initial_fault)


def normalise_name(name: str) -> str:
    """Normalise a name (fault name, realisation name), as a filename or
    YAML key.

    Parameters
    ----------
    name : str
        The name to normalise

    Returns
    -------
    str
        The normalised equivalent of this name. Normalised names are entirely
        lower case, and all non-alphanumeric characters are replaced with "_".
    """
    return re.sub(r"[^A-z0-9]", "_", name.lower())

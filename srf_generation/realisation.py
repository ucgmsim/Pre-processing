#!/usr/bin/env python3
import dataclasses
from pathlib import Path
from typing import Tuple

import yaml
from nshmdb import fault
from nshmdb.fault import Fault

from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)


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


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    dt: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    velocity_model: str
    faults: dict[str, RealisationFault]


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
                planes=[fault.FaultPlane(**params) for params in fault_obj["planes"]],
            )
            for name, fault_obj in faults_object.items()
        }

        for name, fault_obj in faults_object.items():
            if fault_obj["parent"]:
                faults[name].parent = faults[fault_obj["parent"]]

        return Realisation(
            name=raw_yaml_data["name"],
            type=raw_yaml_data["type"],
            dt=raw_yaml_data["dt"],
            genslip_seed=raw_yaml_data["genslip_seed"],
            srfgen_seed=raw_yaml_data["srfgen_seed"],
            genslip_version=raw_yaml_data["genslip_version"],
            faults=faults,
            velocity_model=DEFAULT_1D_VELOCITY_MODEL_PATH,
        )

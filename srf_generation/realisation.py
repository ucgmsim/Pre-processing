#!/usr/bin/env python3
from srf_generation import fault
import dataclasses
from pathlib import Path
import yaml
from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)


@dataclasses.dataclass
class FaultJump:
    parent: fault.Fault
    jump_location_lat: float
    jump_location_lon: float


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    dt: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    velocity_model: str
    faults: dict[str, fault.Fault]


def read_realisation(realisation_filepath: Path) -> Realisation:
    with open(realisation_filepath, "r", encoding="utf-8") as realisation_file:
        raw_yaml_data = yaml.safe_load(realisation_file)
        faults_object = raw_yaml_data["faults"]
        faults = {
            name: fault.Fault(
                name=fault_obj["name"],
                tect_type=fault_obj["tect_type"],
                segments=[
                    fault.FaultSegment(**params) for params in fault_obj["segments"]
                ],
                shyp=fault_obj["shyp"],
                dhyp=fault_obj["dhyp"],
                magnitude=fault_obj["magnitude"],
                parent_jump_coords=(
                    tuple(fault_obj["parent_jump_coords"])
                    if fault_obj["parent_jump_coords"]
                    else None
                ),
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

"""
This module provides classes and functions for working with realisations.

Classes:
    - FaultJump: Represents a jump across two faults.
    - RealisationFault: Represents a fault in a realisation.
    - Realisation: Represents a complete realisation.

Functions:
    - read_realisation: Reads a realisation from a YAML file.
    - topologically_sorted_faults: Generates faults in topologically sorted order.
    - normalise_name: Normalises a name for use in filenames or YAML keys.
"""

import collections
import dataclasses
import re
from pathlib import Path
from typing import Generator, Optional, Tuple

import numpy as np
import yaml

from nshmdb import fault
from nshmdb.fault import Fault
from qcore import coordinates
from qcore.uncertainties import distributions
from srf_generation.source_parameter_generation.common import \
    DEFAULT_1D_VELOCITY_MODEL_PATH


@dataclasses.dataclass
class RealisationFault(Fault):
    """
    Represents a fault in a Realisation.

    Parameters
    ----------
    name : str
        The name of the fault.
    tect_type : str
        The tectonic type of the fault.
    magnitude : float, optional
        The magnitude of the fault.
    parent : RealisationFault, optional
        The parent fault of this fault. Will be None if the fault is the
        initial fault in the rupture.
    parent_jump_coords : tuple, optional
        Coordinates of the parent jump.
    shyp : float, optional
        Distance of hypocentre (along strike).
    dhyp : float, optional
        Distance of hypocenter (along dip).
    """

    magnitude: Optional[float] = None
    parent: Optional["RealisationFault"] = None
    parent_jump_coords: Optional[Tuple[float, float]] = None
    shyp: Optional[float] = None
    dhyp: Optional[float] = None

    def random_fault_coordinates(self) -> Tuple[float, float]:
        """Generate random hypocentre coordinates for the fault.

        Returns
        -------
        tuple
            Random coordinates (shyp, dhyp) for the fault.
        """
        dhyp = self.widths()[0] * distributions.truncated_weibull(1)
        shyp = distributions.rand_shyp() * np.sum(self.lengths())
        return (shyp, dhyp)

    def expected_fault_hypocentre_coordinates(self) -> Tuple[float, float]:
        """Calculate the expected hypocentre coordinates for this fault.

        Returns
        -------
        tuple[float, float]
            Expected coordinates (shyp, dhyp) for the fault.
        """
        dhyp = self.widths()[0] * distributions.truncated_weibull_expected_value(1)
        shyp = 0  # distributions.rand_shyp() is a normal distribution with mean = 0
        return (shyp, float(dhyp))


@dataclasses.dataclass
class Realisation:
    """Represents a complete realisation.

    Parameters
    ----------
    name : str
        The name of the realisation.
    type : int
        The type of the realisation.
    dt : float
        The timestep of the realisation for simulation.
    magnitude : float
        The magnitude of the realisation.
    genslip_seed : int
        The genslip seed of the realisation.
    genslip_version : str
        The genslip version to use for realisation SRF generation.
    srfgen_seed : int
        The srfgen seed of the realisation.
    velocity_model : str
        The velocity model for the realisation.
    faults : dict
        Dictionary of faults associated with the realisation.
    """

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
        """Get the initial fault in the realisation.

        Returns
        -------
        RealisationFault
            The initial fault in the realisation.
        """
        return next(fault for fault in self.faults.values() if not fault.parent)


def read_realisation(realisation_filepath: Path) -> Realisation:
    """Read a realisation from a YAML file.

    Parameters
    ----------
    realisation_filepath : Path
        Path to the YAML file containing realisation data.

    Returns
    -------
    Realisation
        The parsed realisation object.
    """
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
                        coordinates.wgs_depth_to_nztm(np.array(params["corners"])),
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
            magnitude=raw_yaml_data["magnitude"],
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
    """Generate faults in topologically sorted order.

    Parameters
    ----------
    faults : list[RealisationFault]
        List of RealisationFault objects.

    Yields
    ------
    RealisationFault
        The next fault in the topologically sorted order.
    """
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

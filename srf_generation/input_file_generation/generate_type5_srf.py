import dataclasses
from pathlib import Path
from typing import Dict, List

import numpy as np
import pyproj
import yaml

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = Transformer.from_crs(WGS_CODE, NZTM_CODE)
# Convert x, y to lat, lon
NZTM2WGS = Transformer.from_crs(NZTM_CODE, WGS_CODE)


@dataclasses.dataclass
class FaultSegment:
    strike: float
    rake: float
    dip: float
    dtop: float
    dbottom: float
    length: float
    width: float
    dip_dir: float
    clon: float
    clat: float

    def corners(self) -> np.ndarray:
        # assuming that the bottom-left is the origin, the dip is zero and the
        # strike is also zero then the bounds of the plane are easily found.
        corners = np.array(
            [
                [0, 0, 0],  # bottom-left
                [self.width, 0, 0],  # bottom-right
                [self.width, self.length, 0],  # top-right
                [0, self.length, 0],  # top-left
            ]
        ).T
        dip_rad = np.radians(self.dip)
        strike_rad = np.radians(self.strike)

        #
        #   +------+       |    x------x
        #   |      |       |   /      /
        #   |      | --->  |  /      /
        #   |      |       | /      /
        #   |      |       |/      /
        #   +------+       +------+
        #
        #
        #                  |
        #                  |
        #   -------- --->  |
        #                   \
        #                    \
        #                     \

        rot_mat_dip_y = np.array(
            [np.cos(dip_rad), 0, -np.sin(dip_rad)],
            [0, 1, 0],
            [np.sin(dip_rad), 0, np.cos(dip_rad)],
        )
        rot_mat_strike_z = np.array(
            [
                [np.cos(strike_rad), np.sin(strike_rad), 0],
                [-np.sin(strike_rad), np.cos(strike_rad), 0],
                [0, 0, 1],
            ]
        )

        # This now contains the corners relative to the bottom left of fault segment.
        corners_rel_bottom_left = (rot_mat_strike_z @ (rot_mat_dip_y @ corners)).T

        centre_rel_bottom_left_projected_axis_aligned = np.array(
            [self.width * np.cos(dip_rad) / 2, self.length / 2]
        )
        strike_rot = np.array(
            [
                [np.cos(strike_rad), np.sin(strike_rad)],
                [-np.sin(strike_rad), np.cos(strike_rad)],
            ]
        )
        centre_rel_bottom_left_projected = (
            strike_rot @ centre_rel_bottom_left_projected_axis_aligned
        )
        centre_rel_bottom_left_projected = np.append(
            centre_rel_bottom_left_projected, 0
        )
        centre = np.array([*WGS2NZTM.transform(self.clon, self.clat), 0])
        bottom_left = centre - centre_rel_bottom_left_projected
        return corners_rel_bottom_left + bottom_left


Fault = List[FaultSegment]


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    magnitude: float
    moment: float
    dt: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    initial_fault: str
    shypo: float
    dhypo: float
    faults: Dict[str, Fault]


def get_segment_corners(segment: FaultSegment) -> np.ndarray:
    pass


def read_realisation(realisation_filepath: Path) -> Realisation:
    with open(realisation_filepath, "r", encoding="utf-8") as realisation_file:
        raw_yaml_data = yaml.safe_load(realisation_file)
        faults = {
            name: [FaultSegment(**params) for params in subfaults]
            for name, subfaults in raw_yaml_data.pop("faults").items()
        }
        return Realisation(**raw_yaml_data, faults=faults)


if __name__ == "__main__":
    print(read_realisation("/home/jake/src/Pre-processing/test.yaml"))

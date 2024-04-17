import dataclasses
from pathlib import Path
from typing import Dict, List

import numpy as np
import pyproj
import yaml

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
# Convert x, y to lat, lon
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)

KM_CONVERSION_FACTOR = 1000


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

    @property
    def width_km(self):
        return self.width * KM_CONVERSION_FACTOR

    @property
    def length_km(self):
        return self.length * KM_CONVERSION_FACTOR

    @property
    def strike_adjusted(self):
        return self.strike

    def segment_coordinates_to_nztm(
        self, segment_coordinates: np.ndarray
    ) -> np.ndarray:
        """Convert segment coordinates to NZTM global coordinates.

        Parameters
        ----------
        segment_coordinates : np.ndarray
            A (n x 2) matrix of coordinates to convert. Segment coordinates are
            2D coordinates (x, y) given for a fault segment (a plane), where x
            represents displacement along the length of the fault, and y
            displacement along the width of the fault (see diagram below). The
            origin for segment coordinates is the top-left of the fault.

                                     +x
                     0,0 -------------------------->
                        +---------------------------+
                      | |        < length >         |
                      | |                           |  ^
                   +y | |                           | width
                      | |                           |  v
                      v |                           |
                        +---------------------------+
                                                      1,1

        Returns
        -------
        np.ndarray
            An n x 3 matrix of NZTM transformed coordinates.

        Examples
        --------
        FIXME: Add docs.

        """

        # segment coordinates are two-dimensional, so we embed them in three space to begin.

        scaling_factor = np.array([self.width_km, self.length_km])
        segment_coordinates = (segment_coordinates * scaling_factor).T
        segment_coordinates = np.append(
            segment_coordinates,
            np.zeros_like(segment_coordinates[0]).reshape((1, -1)),
            axis=0,
        )
        dip_rad = np.radians(self.dip)
        strike_rad = np.radians(self.strike_adjusted)

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
            [
                [np.cos(dip_rad), 0, -np.sin(dip_rad)],
                [0, 1, 0],
                [np.sin(dip_rad), 0, np.cos(dip_rad)],
            ]
        )
        rot_mat_strike_z = np.array(
            [
                [np.cos(strike_rad), np.sin(strike_rad), 0],
                [-np.sin(strike_rad), np.cos(strike_rad), 0],
                [0, 0, 1],
            ]
        )
        transformation_matrix = rot_mat_strike_z @ rot_mat_dip_y

        centre_rel_bottom_left_projected_axis_aligned = np.array(
            [self.width_km * np.cos(dip_rad) / 2, self.length_km / 2]
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
        centre = np.array([*reversed(WGS2NZTM.transform(self.clat, self.clon)), 0])
        bottom_left = centre - centre_rel_bottom_left_projected
        return (transformation_matrix @ segment_coordinates).T + bottom_left

    def centroid(self) -> np.ndarray:
        """Returns the centre of the fault segment.

        Returns
        -------
        np.ndarray
            A 1 x 3 dimensional vector representing the centroid of the fault
            plane in NZTM coordinates.

        """

        return self.segment_coordinates_to_nztm(np.array([[1 / 2, 1 / 2]]))

    def corners(self) -> np.ndarray:
        """Get the corners of the fault plan

        Returns
        -------
        np.ndarray
            A 4 x 3 dimensional matrix, where each row is a corner of the fault
            plane specified in NZTM coordinates. The corners are returned in
            clockwise orientation starting from the top-left.

        Examples
        --------
        FIXME: Add docs.

        """
        # assuming that the bottom-left is the origin, the dip is zero and the
        # strike is also zero then the bounds of the plane are easily found.
        corners = np.array(
            [
                [0, 0],  # top-left
                [1, 0],  # top-right
                [1, 1],  # bottom-right
                [0, 1],  # bottom-left
            ]
        )
        return self.segment_coordinates_to_nztm(corners)


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

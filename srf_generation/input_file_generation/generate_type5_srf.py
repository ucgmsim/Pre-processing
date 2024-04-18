import dataclasses
from pathlib import Path
from typing import Dict, List

import numpy as np
import pyproj
import yaml
import scipy as sp

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
                      | |        < width >          |
                      | |                           |  ^
                   +y | |                           | length
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
        segment_coordinates = np.append(
            segment_coordinates,
            np.zeros((segment_coordinates.shape[0], 1)),
            axis=1,
        )

        rotation_dip = sp.spatial.transform.Rotation.from_rotvec(
            [0, self.dip, 0], degrees=True
        )
        rotation_strike = sp.spatial.transform.Rotation.from_rotvec(
            [0, 0, self.strike], degrees=True
        )
        centroid_proj = rotation_strike.apply(
            rotation_dip.apply(np.array([self.width_km / 2, self.length_km / 2, 0]))
        )
        centroid_proj[2] = 0
        print(centroid_proj)
        centroid_proj_final = np.array(
            [*reversed(WGS2NZTM.transform(self.clat, self.clon)), 0]
        )
        print(centroid_proj_final)
        centroid_displacement = centroid_proj_final - centroid_proj
        return (
            rotation_strike.apply(rotation_dip.apply(segment_coordinates))
            + centroid_displacement
        )

    def centroid(self) -> np.ndarray:
        """Returns the centre of the fault segment.

        Returns
        -------
        np.ndarray
            A 1 x 3 dimensional vector representing the centroid of the fault
            plane in NZTM coordinates.

        """

        return self.segment_coordinates_to_nztm(
            np.array([[self.width_km / 2, self.length_km / 2]])
        )

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
                [0, self.length_km],  # top-right
                [self.width_km, self.length_km],  # bottom-right
                [self.width_km, 0],  # bottom-left
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

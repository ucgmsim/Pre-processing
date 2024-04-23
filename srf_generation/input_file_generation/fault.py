#!/usr/bin/env python3
import numpy as np
import scipy as sp
import dataclasses
import pyproj
import qcore.geo


# resolution for geometry in width and length, in kilometers

KM_TO_M = 1000
SUBDIVISION_RESOLUTION_KM = 0.1


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
    def width_m(self):
        return self.width * KM_TO_M

    @property
    def length_m(self):
        return self.length * KM_TO_M

    def length_subdivisions(self):
        return int(np.round(self.length / SUBDIVISION_RESOLUTION_KM))

    def width_subdivisions(self):
        return int(np.round(self.width / SUBDIVISION_RESOLUTION_KM))

    def segment_coordinates_to_global_coordinates(
        self, segment_coordinates: np.ndarray
    ) -> np.ndarray:
        """Convert segment coordinates to nztm global coordinates.

        Parameters
        ----------
        segment_coordinates : np.ndarray
            Segment coordinates to convert. Segment coordinates are
            2D coordinates (x, y) given for a fault segment (a plane), where x
            represents displacement along the length of the fault, and y
            displacement along the width of the fault (see diagram below). The
            origin for segment coordinates is the centre of the fault.

                                     +x
               -1/2,-1/2 -------------------------->
                        +---------------------------+
                      | |        < width >          |
                      | |                           |  ^
                   +y | |                           | length
                      | |                           |  v
                      v |                           |
                        +---------------------------+
                                                      1/2,1/2

        Returns
        -------
        np.ndarray
            An 3d-vector of (lat, lon, depth) transformed coordinates.

        Examples
        --------
        FIXME: Add docs.

        """
        projected_width = self.width * np.cos(np.radians(self.dip))
        dip_coord_dir = (
            self.dip_dir if segment_coordinates[0] > 0 else 180 + self.dip_dir
        )
        width_shift_lat, width_shift_lon = qcore.geo.ll_shift(
            self.clat,
            self.clon,
            projected_width * np.abs(segment_coordinates[0]),
            dip_coord_dir,
        )
        strike_coord_dir = (
            self.strike if segment_coordinates[1] > 0 else 180 + self.strike
        )
        length_shift_lat, length_shift_lon = qcore.geo.ll_shift(
            width_shift_lat,
            width_shift_lon,
            self.length * np.abs(segment_coordinates[1]),
            strike_coord_dir,
        )
        depth = (
            (segment_coordinates[0] + 1 / 2)
            * self.width
            * -np.sin(np.radians(self.dip))
        )
        return np.array([length_shift_lat, length_shift_lon, depth])

    def centroid(self) -> np.ndarray:
        """Returns the centre of the fault segment.

        Returns
        -------
        np.ndarray
            A 1 x 3 dimensional vector representing the centroid of the fault
            plane in (lat, lon, depth) format.

        """

        return self.segment_coordinates_to_global_coordinates(np.array([0, 0]))

    def corners(self) -> np.ndarray:
        """Get the corners of the fault plan

        Returns
        -------
        np.ndarray
            A 4 x 3 dimensional matrix, where each row is a corner of the fault
            plane specified in (lat, lon, depth) format. The corners are returned in
            clockwise orientation starting from the top-left.

        Examples
        --------
        FIXME: Add docs.

        """
        # assuming that the bottom-left is the origin, the dip is zero and the
        # strike is also zero then the bounds of the plane are easily found.
        corners = np.array(
            [
                [-1 / 2, -1 / 2],  # top-left
                [-1 / 2, 1 / 2],  # top-right
                [1 / 2, 1 / 2],  # bottom-right
                [1 / 2, -1 / 2],  # bottom-left
            ]
        )
        return np.array(
            [
                self.segment_coordinates_to_global_coordinates(corner)
                for corner in corners
            ]
        )


@dataclasses.dataclass
class Fault:
    name: str
    tect_type: str
    segments: list[FaultSegment]

    def area(self) -> float:
        return sum(segment.width * segment.length for segment in self.segments)

    def widths(self) -> np.ndarray:
        return np.array([seg.width for seg in self.segments])

    def lengths(self) -> np.ndarray:
        return np.array([seg.length for seg in self.segments])

    def number_length_subdivisions(self):
        return self.lengths() / SUBDIVISION_RESOLUTION_KM

    def number_width_subdivisions(self):
        return self.widths() / SUBDIVISION_RESOLUTION_KM

    def corners(self):
        return np.array([segment.corners() for segment in self.segments])

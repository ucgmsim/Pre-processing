#!/usr/bin/env python3
import dataclasses
from typing import Any

import numpy as np
import pyproj
import qcore.geo
import scipy as sp
from qcore.uncertainties import distributions

# resolution for geometry in width and length, in kilometers

KM_TO_M = 1000
SUBDIVISION_RESOLUTION_KM = 0.1

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)


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

    @property
    def projected_width(self):
        return self.length * np.cos(np.radians(self.dip))

    @property
    def projected_width_m(self):
        return self.projected_width * KM_TO_M

    @property
    def depth_m(self):
        return self.length * KM_TO_M

    def length_subdivisions(self):
        return int(np.round(self.length / SUBDIVISION_RESOLUTION_KM))

    def width_subdivisions(self):
        return int(np.round(self.width / SUBDIVISION_RESOLUTION_KM))

    def frame(self):
        strike_direction = np.array(
            [np.cos(np.radians(self.strike)), np.sin(np.radians(self.strike))]
        )
        dip_direction = np.array(
            [
                np.cos(np.radians(self.dip_dir)),
                np.sin(np.radians(self.dip_dir)),
            ]
        )
        return np.array([strike_direction, dip_direction])

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
        depth = (
            (segment_coordinates[0] + 1 / 2)
            * self.width_m
            * -np.sin(np.radians(self.dip))
        )
        dip_coord_dir = (
            self.dip_dir if segment_coordinates[0] > 0 else 180 + self.dip_dir
        )
        projected_width = self.width * np.cos(np.radians(self.dip))
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
        return np.array([length_shift_lat, length_shift_lon, depth])

    def global_coordinates_to_segment_coordinates(
        self,
        global_coordinates: np.ndarray,
    ) -> np.ndarray:
        """Convert coordinates (lat, lon, depth) to segment coordinates (x, y).

        See segment_coordinates_to_global_coordinates for a description of segment
        coordinates.

        Parameters
        ----------
        global_coordinates : np.ndarray
            Global coordinates to convert.

        Returns
        -------
        np.ndarray
            The segment coordinates (x, y) representing the position of
            global_coordinates on the fault segment.

        Examples
        --------
        FIXME: Add docs.

        """
        # Ok how about the stupidest solution ever...

        def f(x):
            return (
                self.segment_coordinates_to_global_coordinates(x[:2])
                - global_coordinates
            )

        return sp.optimize.root(f, np.array([0, 0, 0])).x[:2]

    def coordinate_in_segment(self, global_coordinates: np.ndarray):
        segment_coordinates = self.global_coordinates_to_segment_coordinates(
            global_coordinates
        )
        return np.all(
            np.logical_or(
                np.abs(segment_coordinates) < 1 / 2,
                np.isclose(np.abs(segment_coordinates), 1 / 2, atol=1e-4),
            )
        )

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
        return [
            self.segment_coordinates_to_global_coordinates(corner) for corner in corners
        ]


@dataclasses.dataclass
class Fault:
    name: str
    tect_type: str
    segments: list[FaultSegment]
    shyp: float
    dhyp: float
    parent_jump_coords: (float, float, float) = None
    magnitude: float = None
    parent: Any = None

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

    def hypocentre_wgs_to_fault_coordinates(self, global_coordinates: np.ndarray):
        running_length = 0
        midpoint = np.sum(self.lengths()) / 2
        for segment in self.segments:
            if segment.coordinate_in_segment(global_coordinates):
                segment_coordinates = segment.global_coordinates_to_segment_coordinates(
                    global_coordinates
                )
                strike_length = segment_coordinates[0] + 1 / 2
                dip_length = segment_coordinates[1] + 1 / 2
                return (
                    running_length + strike_length * segment.length - midpoint,
                    max(dip_length * segment.width, 0),
                )
            running_length += segment.length
        raise ValueError("Specified coordinates not contained on fault.")

    def fault_coordinates_to_wgsdepth_coordinates(
        self, fault_coordinates: np.ndarray
    ) -> np.ndarray:
        midpoint = np.sum(self.lengths()) / 2
        remaining_length = fault_coordinates[0] + midpoint
        for segment in self.segments:
            if remaining_length < segment.length:
                return segment.segment_coordinates_to_global_coordinates(
                    np.array([remaining_length / segment.length - 1 / 2,
                              fault_coordinates[1] / segment.width - 1 / 2]),
                )
            remaining_length -= segment.length
        raise ValueError("Specified fault coordinates out of bounds.")

    def random_fault_coordinates(self) -> (float, float):
        weibull_scale = 0.612
        dhyp = self.widths()[0] * sp.stats.weibull_min(
            3.353,
            0,
            1,
            scale=weibull_scale
        ).rvs(1)[0]
        shyp = distributions.rand_shyp() * np.sum(self.lengths())
        return (shyp, dhyp)

    def expected_fault_coordinates(self) -> (float, float):
        weibull_scale = 0.612
        dhyp = self.widths()[0] * sp.stats.weibull_min(
            3.353,
            0,
            1,
            scale=weibull_scale
        ).expect()
        shyp = 0
        return (shyp, dhyp)

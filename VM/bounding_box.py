import dataclasses

import numpy as np
from qcore import geo


@dataclasses.dataclass
class BoundingBox:
    """BoundingBox."""

    corners: np.ndarray

    @property
    def origin(self):
        """origin."""
        return np.mean(self.corners, axis=0)

    @property
    def extent_x(self):
        """extent_x."""
        return np.linalg.norm(np.linalg.norm(self.corners[2] - self.corners[1]) / 1000)

    @property
    def extent_y(self):
        """extent_y."""
        return np.linalg.norm(np.linalg.norm(self.corners[1] - self.corners[0]) / 1000)

    @property
    def bearing(self):
        """bearing."""

        north_direction = np.array([1, 0, 0])
        up_direction = np.array([0, 0, 1])
        horizontal_direction = np.append(self.corners[1] - self.corners[0], 0)
        return geo.oriented_bearing_wrt_normal(
            north_direction, horizontal_direction, up_direction
        )

    @property
    def area(self):
        """Return the box area (km^2)"""
        return self.extent_x * self.extent_y

    @property
    def polygon(self):
        """Return a shapely geometry for the bounding box."""
        return shapely.Polygon(
            np.append(self.corners, np.atleast_2d(self.corners[0]), axis=0)
        )


def axis_aligned_bounding_box(points: np.ndarray) -> BoundingBox:
    """Return an axis-aligned bounding box containing points.

    Parameters
    ----------
    points : np.ndarray
        The points to bound.

    Returns
    -------
    np.ndarray
        A numpy array of shape (4, 2) containing defining the corners of
        the bounding box. The order of points is bottom-left, bottom-right,
        top-right, top-left.
    """
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    corners = np.array([[min_x, min_y], [max_x, min_y], [max_x, max_y], [min_x, max_y]])
    return BoundingBox(corners)


def rotation_matrix(angle: float) -> np.ndarray:
    """Return the 2D rotation matrix for a given angle.

    Parameters
    ----------
    angle : float
        The angle to rotate by in radians.

    Returns
    -------
    np.ndarray
        The 2x2 matrix rotating points by angle radians counter-clockwise.
    """
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


def minimum_area_bounding_box(points: np.ndarray) -> BoundingBox:
    """Return the smallest rectangle bounding points. The rectangle may
    be rotated.

    Parameters
    ----------
    points : np.ndarray
        The points to bound.

    Returns
    -------
    BoundingBox
        The bounding box.
    """
    convex_hull = sp.spatial.ConvexHull(points).points
    segments = np.array(
        [
            convex_hull[(i + 1) % len(convex_hull)] - convex_hull[i]
            for i in range(len(convex_hull))
        ]
    )
    rotation_angles = -np.arctan2(segments[:, 1], segments[:, 0])

    bounding_boxes = [
        axis_aligned_bounding_box(convex_hull @ rotation_matrix(angle).T)
        for angle in rotation_angles
    ]

    rotation_angle, minimum_bounding_box = min(
        zip(rotation_angles, bounding_boxes), key=lambda rot_box: rot_box[1].area
    )
    return BoundingBox(
        minimum_bounding_box.corners @ rotation_matrix(-rotation_angle).T
    )

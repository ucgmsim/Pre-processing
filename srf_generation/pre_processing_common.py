import numpy as np

from srf_generation.source_parameter_generation.uncertainties.common import ONE_DEG_LAT


def calculate_corners(dip, x_pos, y_pos, lat, lon, strike):
    # use cartesian coordinate system to define the along strike and downdip
    # locations taking the center of the fault plane as (x,y)=(0,0)

    nx = len(x_pos)
    ny = len(y_pos)
    # now use a coordinate transformation to go from fault plane to North and
    # East cartesian plane (again with (0,0) ==(0,0)  )
    y_pos_surf_proj = y_pos * np.cos(np.radians(dip))
    azimuth = np.radians(strike - 90)
    rot_matrix = np.array(
        [[np.cos(azimuth), np.sin(azimuth)], [-np.sin(azimuth), np.cos(azimuth)]]
    )
    east_loc_relative, north_loc_relative = np.dot(
        rot_matrix, [np.tile(x_pos, ny), y_pos_surf_proj.repeat(ny)]
    ).reshape(2, nx, ny)
    # now use a coordinate transformation to go from the North East cartesian
    # plane to the spherical earth WGS84 coordinate system
    lats = lat + north_loc_relative / ONE_DEG_LAT
    lons = lon + (east_loc_relative / ONE_DEG_LAT) * 1 / np.cos(np.radians(lat))
    return lats, lons


def get_hypocentre(lat, lon, shypo, dhypo, strike, dip):
    """
    Same logic as corners, for a single point.
    """
    azimuth = np.radians(strike - 90)
    rot_matrix = np.array(
        [[np.cos(azimuth), np.sin(azimuth)], [-np.sin(azimuth), np.cos(azimuth)]]
    )

    y_pos_surf_proj_hyp = -dhypo * np.cos(np.radians(dip))
    A = np.dot(rot_matrix, [[shypo], [y_pos_surf_proj_hyp]])

    lat_hyp = lat + A[1][0] / ONE_DEG_LAT
    lon_hyp = lon + (A[0][0] / ONE_DEG_LAT) / np.cos(np.radians(lat))

    return lon_hyp, lat_hyp

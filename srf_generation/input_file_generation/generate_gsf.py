import json
from typing import TextIO

import numpy as np
import pyproj
from srf_generation import fault

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)


def wgsdepth_to_nztm(wgsdepthcoordinates: np.ndarray) -> np.ndarray:
    nztm_coords = np.array(
        WGS2NZTM.transform(wgsdepthcoordinates[:, 0], wgsdepthcoordinates[:, 1]),
    ).T
    return np.append(nztm_coords, wgsdepthcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def nztm_to_wgsdepth(nztmcoordinates: np.ndarray) -> np.ndarray:
    wgs_coords = np.array(
        NZTM2WGS.transform(nztmcoordinates[:, 0], nztmcoordinates[:, 1]),
    ).T
    return np.append(wgs_coords, nztmcoordinates[:, 2].reshape((-1, 1)), axis=-1)


def gridpoints_along_direction(length: float, resolution: float) -> int:
    return int(np.round(length / resolution + 2))


def fault_segment_to_meshgrid(
    segment: fault.FaultSegment, resolution=100
) -> np.ndarray:
    corners = wgsdepth_to_nztm(np.array(segment.corners()))
    origin = corners[0]
    x_upper = corners[1]
    length_x = segment.length_m
    y_bottom = corners[-1]
    length_y = segment.width_m
    nx = gridpoints_along_direction(length_x, resolution)
    ny = gridpoints_along_direction(length_y, resolution)
    x, x_step = np.linspace(0, length_x, nx, retstep=True)
    y, y_step = np.linspace(0, length_y, ny, retstep=True)
    xv, yv = np.meshgrid(x, y)
    coordinates = np.vstack([xv.ravel(), yv.ravel()])
    transformation_matrix = np.vstack(
        [(x_upper - origin) / length_x, (y_bottom - origin) / length_y]
    ).T
    nztm_meshgrid = (transformation_matrix @ coordinates).T
    nztm_meshgrid += origin
    return nztm_to_wgsdepth(nztm_meshgrid)


def write_fault_to_gsf_file(
    gsf_file_handle: TextIO, fault: fault.Fault, resolution=100
):
    meshgrids = [
        fault_segment_to_meshgrid(segment, resolution) for segment in fault.segments
    ]
    number_of_points = sum(meshgrid.shape[0] for meshgrid in meshgrids)
    number_of_strike_gridpoints = sum(gridpoints_along_direction(segment.length_m, resolution) for segment in fault.segments)
    number_of_dip_gridpoints  = gridpoints_along_direction(fault.segments[0].width_m, resolution)
    fault_length = np.sum(fault.lengths())
    fault_width = fault.segments[0].width
    gsf_file_handle.write(f'# nstk= {number_of_strike_gridpoints} ndip= {number_of_dip_gridpoints}\n')
    gsf_file_handle.write(f'# flen= {fault_length:10.4f} fwid={fault_width:10.4f}\n')
    gsf_file_handle.write("# LON  LAT  DEP(km)  SUB_DX  SUB_DY  LOC_STK  LOC_DIP  LOC_RAKE  SLIP(cm)  INIT_TIME  SEG_NO\n")
    gsf_file_handle.write(f"{number_of_points}\n")
    for i, meshgrid in enumerate(meshgrids):
        segment = fault.segments[i]
        strike_step = segment.length / gridpoints_along_direction(
            segment.length_m, resolution
        )
        dip_step = segment.width / gridpoints_along_direction(
            segment.width_m, resolution
        )
        for point in meshgrid:
            gsf_file_handle.write(
                f"{point[1]:11.5f} {point[0]:11.5f} {point[2] / 1000:11.5e} {strike_step:11.5e} {dip_step:11.5e} {segment.strike:6.1f} {segment.dip:6.1f} {segment.rake:6.1f} {-1.0:8.2f} {-1.0:8.2f} {i:3d}\n"
            )

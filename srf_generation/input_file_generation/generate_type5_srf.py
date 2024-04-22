import dataclasses
from pathlib import Path
import matplotlib.pyplot as plt
import multiprocessing
import rupture_propogation
import qcore.geo

from srf_generation.source_parameter_generation.common import (
    DEFAULT_1D_VELOCITY_MODEL_PATH,
)
import numpy as np
import pyproj
import scipy as sp
import yaml
from qcore import binary_version
import subprocess
import tempfile

WGS_CODE = 4326
NZTM_CODE = 2193
# Convert lat, lon to x, y
WGS2NZTM = pyproj.Transformer.from_crs(WGS_CODE, NZTM_CODE)
UTM_CRS = pyproj.CRS(proj="utm", zone=11, ellps="WGS84")
WGS2UTM = pyproj.Transformer.from_crs(WGS_CODE, UTM_CRS)
TRANSFORMER_MAP = {"nztm": WGS2NZTM, "utm": WGS2UTM}
# Convert x, y to lat, lon
NZTM2WGS = pyproj.Transformer.from_crs(NZTM_CODE, WGS_CODE)
# resolution for geometry in width and length, in kilometers
SUBDIVISION_RESOLUTION_KM = 0.1

KM_TO_M = 1000
FAULTSEG2GSFDIPDIR = "fault_seg2gsf_dipdir"


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
        return self.width * KM_TO_M

    @property
    def length_km(self):
        return self.length * KM_TO_M

    @property
    def strike_adjusted(self):
        return self.strike

    def length_subdivisions(self):
        return int(np.round(self.length / SUBDIVISION_RESOLUTION_KM))

    def width_subdivisions(self):
        return int(np.round(self.width / SUBDIVISION_RESOLUTION_KM))

    def segment_coordinates_to_global_coordinates(
        self, segment_coordinates: np.ndarray, coordinate_system="nztm"
    ) -> np.ndarray:
        """Convert segment coordinates to nztm global coordinates.

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
            An n x 3 matrix of nztm transformed coordinates.

        Examples
        --------
        FIXME: Add docs.

        """
        # The approach to transforming these coordinates is different to the
        # approach taken in past. We are going to start with the simplest
        # possible case of a flat fault and then transform the coordinates to
        # match the fault we
        # actually want the coordinates of.

        # Suppose that we fault with a strike and dip of exactly 0, a length l
        # and a width w. Then the segment coordinates (x, y) correspond to the
        # 3D coordinates (x, y, 0). So here we just add that zero coordinate on.
        segment_coordinates = np.append(
            segment_coordinates,
            np.zeros((segment_coordinates.shape[0], 1)),
            axis=1,
        )

        # To transform this flat fault into an arbitrary fault, we simply rotate
        # first through dip and then strike. The dip is always to the right of
        # the strike (i.e. clockwise from the y-axis). Rotation vectors by
        # default rotate anti-clockwise so we need to rotate by -dip.
        rotation_dip = sp.spatial.transform.Rotation.from_rotvec(
            [0, -self.dip, 0], degrees=True
        )
        # Similarly, strike is a clockwise oriented rotation from the z-axis, so
        # we need to rotate by -strike.
        rotation_strike = sp.spatial.transform.Rotation.from_rotvec(
            [0, 0, -self.strike], degrees=True
        )
        # We need to translate the our coordinates to match the position of the
        # fault. The fault definition gives us the position of the projection of
        # the centroid onto the earth's surface in lat-lon coordinates. We will
        # translate the projected centroid to match the centroid given in the
        # fault definition.
        centroid_proj = rotation_strike.apply(
            rotation_dip.apply(np.array([self.width_km / 2, self.length_km / 2, 0]))
        )
        centroid_proj[2] = 0
        # Now we lookup the projected centroid location from the fault
        # definition and construct the displacement vector.
        coordinate_transformer = TRANSFORMER_MAP[coordinate_system]
        centroid_proj_final = np.array(
            [*reversed(coordinate_transformer.transform(self.clat, self.clon)), 0]
        )
        centroid_displacement = centroid_proj_final - centroid_proj
        # So the transformed coordinates is
        # rotation * coordinates + centroid displacement.
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

        return self.segment_coordinates_to_global_coordinates(
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
        return self.segment_coordinates_to_global_coordinates(corners)


@dataclasses.dataclass
class Fault:
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
        return np.array([segment.corners()] for segment in self.segments)


@dataclasses.dataclass
class Realisation:
    name: str
    type: int
    magnitude: float
    # moment: float
    dt: float
    genslip_seed: int
    genslip_version: str
    srfgen_seed: int
    initial_fault: str
    shypo: float
    dhypo: float
    velocity_model: str
    faults: dict[str, Fault]

    def rupture_area(self) -> float:
        return sum(fault.area() for fault in self.faults.values())

    def fault_magnitude_by_area(self, fault: str) -> float:
        return self.magnitude * (self.faults[fault].area() / self.rupture_area())


def read_realisation(realisation_filepath: Path) -> Realisation:
    with open(realisation_filepath, "r", encoding="utf-8") as realisation_file:
        raw_yaml_data = yaml.safe_load(realisation_file)
        faults = {
            name: Fault(
                tect_type=fault["tect_type"],
                segments=[FaultSegment(**params) for params in fault["segments"]],
            )
            for name, fault in raw_yaml_data.pop("faults").items()
        }
        return Realisation(
            **raw_yaml_data,
            faults=faults,
            velocity_model=DEFAULT_1D_VELOCITY_MODEL_PATH,
        )


def type4_fault_gsf(realisation: Realisation, fault: Fault):
    fault_seg_bin = binary_version.get_unversioned_bin(FAULTSEG2GSFDIPDIR)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as input_file, tempfile.NamedTemporaryFile(
        mode="w", delete=False
    ) as gsf_output_file:
        input_file.write(f"{len(fault.segments)}\n")
        for segment in fault.segments:
            input_file.write(
                f"{segment.clon:6f} {segment.clat:6f} {segment.dtop:6f} {segment.strike:.4f} {segment.dip:.4f} {segment.rake:.4f} {segment.length:.4f} {segment.width:.4f} {segment.length_subdivisions():d} {segment.width_subdivisions():d}\n"
            )
        input_file.flush()
        fault_seg_command = [
            fault_seg_bin,
            "read_slip_vals=0",
            f"infile={input_file.name}",
            f"outfile={gsf_output_file.name}",
            f"dipdir={segment.dip_dir}",
        ]
        print(" ".join(fault_seg_command))
        subprocess.run(fault_seg_command, check=True)

    return gsf_output_file.name


def generate_type4_fault_srf(
    realisation: Realisation,
    fault_name: str,
    output_directory: Path,
    hypocentre: np.ndarray,
):
    fault = realisation.faults[fault_name]

    gsf_file_path = type4_fault_gsf(realisation, fault)

    genslip_bin = binary_version.get_genslip_bin(realisation.genslip_version)
    magnitude = realisation.fault_magnitude_by_area(fault_name)

    lengths = fault.lengths()
    nx = int(np.sum(np.round(lengths / SUBDIVISION_RESOLUTION_KM)))
    ny = int(np.round(fault.widths()[0] / SUBDIVISION_RESOLUTION_KM))

    genslip_cmd = [
        genslip_bin,
        "read_erf=0",
        "write_srf=1",
        "read_gsf=1",
        "write_gsf=0",
        f"infile={gsf_file_path}",
        f"mag={magnitude}",
        f"nstk={nx}",
        f"ndip={ny}",
        "ns=1",
        "nh=1",
        f"seed={realisation.genslip_seed}",
        f"velfile={realisation.velocity_model}",
        f"shypo={hypocentre[0]}",
        f"dhypo={hypocentre[1]}",
        f"dt={realisation.dt}",
        "plane_header=1",
        "srf_version=1.0",
        "seg_delay={0}",
        "rvfac_seg=-1",
        "gwid=-1",
        "side_taper=0.02",
        "bot_taper=0.02",
        "top_taper=0.0",
        "rup_delay=0",
        "alpha_rough=0.0",
    ]

    if fault.tect_type == "SUBDUCTION_INTERFACE":
        genslip_cmd.extend(
            [
                "kmodel=-1",
                "xmag_exp=0.5",
                "ymag_exp=0.5",
                "kx_corner=2.5482",
                "ky_corner=2.3882",
                "tsfac_slope=-0.5",
                "tsfac_bzero=-0.1",
                "risetime_coef=1.95",
            ]
        )

    srf_file_path = output_directory / (fault_name + ".srf")
    with open(srf_file_path, "w", encoding="utf-8") as srf_file_handle:
        subprocess.run(
            genslip_cmd, stdout=srf_file_handle, stderr=subprocess.PIPE, check=True
        )

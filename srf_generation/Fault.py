from abc import ABC
from typing import Union

import numpy as np
from qcore import geo
from qcore.nhm import NHMFault

from srf_generation.source_parameter_generation.uncertainties.mag_scaling import (
    MagnitudeScalingRelations,
    get_area,
    get_length,
    get_width,
    mag2mom,
    round_subfault_size,
    lw_2_mw_scaling_relation
)


def fault_factory(fault_type: int):
    return [Type1, Type2, Type3, Type4][fault_type - 1]


class Fault(ABC):
    subfault_spacing = 0.1
    type = 0
    _mag = None
    _moment = None
    name = None
    _latitude = None
    _longitude = None
    _strike = None
    _rake = None
    _dip = None
    _depth = None
    _width = None
    _length = None
    _shypo = None
    _dhypo = None

    @property
    def pid(self):
        return self.name

    @pid.setter
    def pid(self, value):
        self.name = value

    @property
    def mom(self):
        if self._moment is None:
            self._moment = mag2mom(self._mag)
        return self._moment

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        self._latitude = value

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        self._longitude = value

    @property
    def magnitude(self):
        return self._mag

    @magnitude.setter
    def magnitude(self, mag):
        if mag > 11:
            raise ValueError(f"Given mag {mag} is greater than theoretically possible 11")
        self._mag = mag

    @property
    def length(self):
        return self._length

    @property
    def width(self):
        return self._width

    def to_dict(self):
        return {
            "type": self.type,
            "magnitude": self._mag,
            "moment": self.mom,
            "name": self.name,
            "longitude": self.longitude,
            "latitude": self.latitude,
            "strike": self._strike,
            "rake": self._rake,
            "dip": self._dip,
            "depth": self._depth,
        }

    @property
    def rake(self):
        return self._rake

    @rake.setter
    def rake(self, value):
        value = ((value + 180) % 360) - 180
        self._set_rake(value)

    def _set_rake(self, value):
        self._rake = value


class SinglePlaneFault(Fault):
    def __init__(self, name, magnitude, strike, rake, dip):
        self.name = name
        self._mag = magnitude
        self._strike = strike
        self._rake = rake
        self._dip = dip

    _length: float = None
    _width: float = None
    _dtop: float = None
    _dbottom: float = None

    @property
    def strike(self):
        return self._strike

    @strike.setter
    def strike(self, value):
        self._strike = value % 360

    @property
    def dip(self):
        return self._dip

    @dip.setter
    def dip(self, value):
        if 90 < value < 180:
            self.strike = self.strike + 180
            value = 180 - value
        elif value < 0 or value > 180:
            raise ValueError("Invalid dip value")

        self._dip = value

    @property
    def dhypo(self):
        return self.width / 2

    @property
    def dbottom(self):
        return self._dbottom

    @property
    def dtop(self):
        return self._dtop

    @property
    def hypocentre_lonlat(self):
        return self._longitude, self._latitude


class MultiPlaneFault(Fault):
    pass


class PointSourceFault(SinglePlaneFault):

    type = 1

    # vs = 3.2
    # rh0 = 2.44
    # risetime = 0.5
    # stype = "cos"
    # inittime = 0.0
    vs = None
    rh0 = None
    risetime = None
    stype = None
    inittime = None

    def __init__(self, name, latitude, longitude, magnitude, strike, rake, dip, depth):
        self._depth = depth
        self._latitude = latitude
        self._longitude = longitude
        super().__init__(name, magnitude, strike, rake, dip)

    @property
    def dbottom(self):
        return self._depth

    def to_dict(self):
        base_dict = super().to_dict()
        base_dict.update({"type": 1})
        if self.vs is not None:
            base_dict["vs"] = self.vs
        if self.rh0 is not None:
            base_dict["rho"] = self.rh0
        if self.risetime is not None:
            base_dict["risetime"] = self.risetime
        if self.stype is not None:
            base_dict["stype"] = self.stype
        if self.inittime is not None:
            base_dict["inittime"] = self.inittime
        return base_dict


Type1 = PointSourceFault


class FiniteFault(SinglePlaneFault):

    dlen = 0.1
    dwid = 0.1
    shypo = 0.00
    depths = None

    def to_dict(self):
        base_dict: dict = super().to_dict()
        base_dict.update(
            {
                "flen": self.length,
                "dlen": self.dlen,
                "fwid": self.width,
                "dwid": self.dwid,
                "dtop": self.dtop,
                "shypo": self.shypo,
                "dhypo": self.dhypo,
            }
        )
        return base_dict


class Type2(FiniteFault):

    type = 2

    _ratio_override: Union[float, None] = None
    _magnitude_scaling_relation: MagnitudeScalingRelations = None

    def __init__(self, name, latitude, longitude, magnitude, strike, rake, dip, depth):
        self._depth = depth
        self._latitude = latitude
        self._longitude = longitude
        super().__init__(name, magnitude, strike, rake, dip)

    def to_dict(self):
        base_dict: dict = super().to_dict()
        base_dict.update(
            {
                "type": 2,
                "longitude": self._lon_hyp,
                "latitude": self._lat_hyp,
                "mwsr": self._magnitude_scaling_relation.name,
            }
        )
        return base_dict

    def _set_rake(self, value):
        self._rake = value
        if self._magnitude_scaling_relation is not None:
            self._calculate_dimensions()

    def set_ratio_override(self, ratio):
        self._ratio_override = ratio

    def reset_ratio_override(self):
        self._ratio_override = None

    @property
    def magnitude_scaling_relation(self):
        return self._magnitude_scaling_relation

    @magnitude_scaling_relation.setter
    def magnitude_scaling_relation(self, value: MagnitudeScalingRelations):
        self._magnitude_scaling_relation = value
        self._calculate_dimensions()

    @property
    def length(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")

        return self._length

    @property
    def width(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")

        return self._width

    @property
    def dhypo(self):
        return self.width / 2

    @property
    def dbottom(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        return self._dbottom

    @property
    def dtop(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        return self._dtop

    @property
    def hypocentre_lonlat(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        return self._lon_hyp, self._lat_hyp

    def _calculate_dimensions(self):
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")

        area = get_area(self)
        length = get_length(self)
        width = get_width(self)

        if self._ratio_override is not None:
            r = self._ratio_override
        else:
            r = max(length / width, 1)

        self._width = np.sqrt(area / r)
        self._length = r * np.sqrt(area / r)

        self._dtop = self._depth - np.sin(np.radians(self._dip)) * self.width / 2
        shift = min(self._dtop, 0)
        if shift != 0:
            self._dtop = 0
        self._dbottom = (
            self._depth + np.sin(np.radians(self._dip)) * self.width / 2 + shift
        )
        if self._dbottom > 12:
            self._dbottom += 3

        self.ny = int(round(self.width / self.dwid))
        self.nx = int(round(self.length / self.dlen))
        self.dlen = self.length / float(self.nx)
        self.dwid = self.width / float(self.ny)

        # self._calculate_finite_fault_properties()
        # self._calculate_corners()
        self._get_hypocentre(0, -self.dhypo)

    def _get_hypocentre(self, shypo, dhypo):
        """
        Same logic as corners, for a single point.
        """
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        ONE_DEG_LAT = np.radians(6371.0072)
        azimuth = np.radians(self.strike - 90)
        rot_matrix = np.array(
            [[np.cos(azimuth), np.sin(azimuth)], [-np.sin(azimuth), np.cos(azimuth)]]
        )

        y_pos_surf_proj_hyp = -dhypo * np.cos(np.radians(self._dip))
        hypocentre_points = np.dot(rot_matrix, [[shypo], [y_pos_surf_proj_hyp]])

        self._lat_hyp = self.latitude + hypocentre_points[1][0] / ONE_DEG_LAT
        self._lon_hyp = self.longitude + (
            hypocentre_points[0][0] / ONE_DEG_LAT
        ) / np.cos(np.radians(self.latitude))

    def _calculate_corners(self):
        # use cartesian coordinate system to define the along strike and downdip
        # locations taking the center of the fault plane as (x,y)=(0,0)
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        ONE_DEG_LAT = np.radians(6371.0072)

        x_pos = np.arange(self.dlen / 2.0, self.length, self.dlen) - self.length / 2.0
        y_pos = (
            np.arange(self.dwid / 2.0, self.width, self.dwid)[::-1] - self.width / 2.0
        )

        # now use a coordinate transformation to go from fault plane to North and
        # East cartesian plane (again with (0,0) ==(0,0)  )
        y_pos_surf_proj = y_pos * np.cos(np.radians(self._dip))
        azimuth = np.radians(self.strike - 90)
        rot_matrix = np.array(
            [[np.cos(azimuth), np.sin(azimuth)], [-np.sin(azimuth), np.cos(azimuth)]]
        )
        east_loc_relative, north_loc_relative = np.dot(
            rot_matrix, [np.tile(x_pos, self.ny), y_pos_surf_proj.repeat(self.nx)]
        ).reshape((2, self.nx, self.ny))
        # now use a coordinate transformation to go from the North East cartesian
        # plane to the spherical earth WGS84 coordinate system
        self.lats = self.latitude + north_loc_relative / ONE_DEG_LAT
        self.lons = self.longitude + (east_loc_relative / ONE_DEG_LAT) * 1 / np.cos(
            np.radians(self.latitude)
        )

    def _calculate_finite_fault_properties(self):
        """
        Purpose: To create a finite fault geometry based on centroid moment tensor
        solution information and magnitude scaling relationships.
        The use of such finite fault geometry is for first order computation of
        source-to-site distances.

        revisions
        v2 - reoriented fault plane to be consistent with along strike direction
        v3 - added 'argout_emod3d' output and associated commands to output details
             which are input in the finite fault model generation for EMOD3D graves methodology

        assumptions:
        1) The centroid moment tensor represents the centroid of the fault plane,
        located half way along strike and downdip
        2) The fault area is computed from the median of a Moment scaling
        relationship and thus is assumed to be a 'typical' stress drop event for
        the Mw-relationship used.
        Should only be used for small events where the downdip width is smaller
        than the seismogenic depth (i.e. so the assumption of a square fault plane
        is reasonable); depth is reset to zero if above ground values determined
        3) srike in deg measured clockwise from N; rake in deg measured
        anti-clockwise from the strike direction (and in the range
        -180<lambda<180; dip in deg measured downward from horizontal

        Input variables:
        Lat    - the latitude of the focal mechanism (-ve below equator)
        Lon    - the lon of the focal mech (-180<lon<180)
        Depth  - the depth of the focal mech (in km)
        Mw     - the moment magnitude from the Mw tensor soln
        strike - the strike of the focal mech (only 1)
        rake   - rake (not used, but carried forward)
        dip    - the dip angle of the fault plane

        output variables:
        lat       - the latitude of the subfault
        lon       - the lon of the subfault
        depth     - depth of the subfault
        lonAsList - as a 1D array
        """
        if self.magnitude_scaling_relation is None:
            raise ValueError("No magnitude scaling relation given")
        # rounded subfault spacing

        # use cartesian coordinate system to define the along strike and downdip
        # locations taking the center of the fault plane as (x,y)=(0,0)
        y_pos: np.ndarray = np.arange(self.dwid / 2.0, self.width, self.dwid)[
            ::-1
        ] - self.width / 2.0

        depth_loc_relative = (
            (-y_pos * np.sin(np.radians(self._dip)))
            .repeat(self.nx)
            .reshape((self.nx, self.ny))
        )
        shift = np.min(self.dtop + depth_loc_relative, 0)
        self.depths = self.dtop + depth_loc_relative - shift
        if np.any(0 > self.depths):
            raise ValueError(
                "Some points are above the ground. This represents a logic problem and should be referred to the "
                "developers."
            )


class Type3(FiniteFault):
    def __init__(
        self, name, magnitude, lon1, lat1, lon2, lat2, dip, rake, dtop, dbottom, dip_dir
    ):
        self._trace = ((lon1, lat1), (lon2, lat2))
        self._dtop = dtop
        self._dbottom = dbottom
        self._dip_dir = dip_dir

        self._length = round_subfault_size(geo.ll_dist(lon1, lat1, lon2, lat2), magnitude)
        if dbottom >= 12:
            raw_fwid = (dbottom - dtop + 3) / np.sin(np.radians(dip))
        else:
            raw_fwid = (dbottom - dtop) / np.sin(np.radians(dip))
        self._width = round_subfault_size(raw_fwid, magnitude)

        strike = geo.ll_bearing(lon1, lat1, lon2, lat2)

        super().__init__(name, magnitude, strike, rake, dip)

    def to_dict(self):
        base_dict = {
            "clon": geo.ll_mid(*self._trace[0], *self._trace[1])[0],
            "clat": geo.ll_mid(*self._trace[0], *self._trace[1])[1],
            "type": self.type,
            "magnitude": self._mag,
            "moment": self.mom,
            "name": self.name,
            "strike": self._strike,
            "rake": self._rake,
            "dip": self._dip,
            "dtop": self._dtop,
            "dbottom": self._dbottom,
            "length": self._length,
            "width": self._width,
            "dip_dir": self._dip_dir,
        }
        return base_dict


class Type4(MultiPlaneFault):

    type = 4

    def __init__(self, nhm_data: NHMFault):
        self.name = nhm_data.name
        self._dbottom = nhm_data.dbottom
        self.fault_type = nhm_data.fault_type
        self.tectonic_type = nhm_data.tectonic_type
        self._length = nhm_data.length
        self._dip = nhm_data.dip
        self._rake = nhm_data.rake
        self._dtop = nhm_data.dtop
        self._slip_rate = nhm_data.slip_rate
        self._dip_dir = nhm_data.dip_dir

        self._n_planes = len(nhm_data.trace) - 1

        if nhm_data.tectonic_type == "SUBDUCTION_INTERFACE":
            self.mwsr = MagnitudeScalingRelations.SKARLATOUDIS2016

        else:
            self.mwsr = MagnitudeScalingRelations.LEONARD2014

        dummy_plane = Type3(
            nhm_data.name,
            nhm_data.mw,
            *nhm_data.trace[0],
            *nhm_data.trace[1],
            nhm_data.dip,
            nhm_data.rake,
            nhm_data.dtop,
            nhm_data.dbottom,
            nhm_data.dip_dir,
        )
        self._mag = lw_2_mw_scaling_relation(
            nhm_data.length,
            dummy_plane.width,
            self.mwsr,
            nhm_data.rake
        )

        self._planes = []
        for i in range(self._n_planes):
            self._planes.append(
                Type3(
                    nhm_data.name,
                    self._mag,
                    *nhm_data.trace[i],
                    *nhm_data.trace[i + 1],
                    nhm_data.dip,
                    nhm_data.rake,
                    nhm_data.dtop,
                    nhm_data.dbottom,
                    nhm_data.dip_dir,
                )
            )

    @property
    def length(self):
        return sum([f.length for f in self._planes])

    @property
    def width(self):
        return self._planes[0].width

    @property
    def dhypo(self):
        return self._dhypo

    @dhypo.setter
    def dhypo(self, dhypo):
        if dhypo > self.width:
            raise ValueError(
                f"Cannot place hypocentre outside fault plane. dhpyo: {dhypo}, fault width: {self.width}"
            )
        self._dhypo = dhypo
        for sub_plane in self._planes:
            sub_plane._dhypo = dhypo

    @property
    def shypo(self):
        return self._shypo

    @shypo.setter
    def shypo(self, shypo):
        if shypo > self.length:
            raise ValueError(
                f"Cannot place hypocentre outside fault plane. shpyo: {shypo}, fault length: {self.length}"
            )
        self._shypo = shypo
        for sub_plane in self._planes:
            sub_plane._shypo = shypo

    def to_dict(self):

        base_dict = {
            "type": self.type,
            "magnitude": self._mag,
            "moment": self.mom,
            "name": self.name,
            "fault_type": self.fault_type,
            "tect_type": self.tectonic_type,
            "rake": self._rake,
            "dip": self._dip,
            "dtop": self._dtop,
            "dbottom": self._dbottom,
            "length": self._length,
            "plane_count": self._n_planes,
            "slip_rate": self._slip_rate,
            "dip_dir": self._dip_dir,
            "shypo": self.shypo,
            "dhypo": self.dhypo,
        }
        for i, sub_fault in enumerate(self._planes):
            for key, item in sub_fault.to_dict().items():
                base_dict[f"{key}_subfault_{i}"] = item
        return base_dict

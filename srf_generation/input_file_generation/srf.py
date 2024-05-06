#!/usr/bin/env python3

import dataclasses
from pathlib import Path
from typing import TextIO
import re


@dataclasses.dataclass
class SrfSegment:
    elon: float
    elat: float
    nstk: int
    ndip: int
    len: float
    wid: float
    stk: float
    dip: float
    dtop: float
    shyp: float
    dhyp: float


@dataclasses.dataclass
class SrfPoint:
    lon: float
    lat: float
    dep: float
    stk: float
    dip: float
    area: float
    tinit: float
    dt: float
    rake: float
    slip1: float
    slip2: float
    slip3: float
    sr1: list[float]
    sr2: list[float]
    sr3: list[float]


@dataclasses.dataclass
class SrfFile:
    header: list[SrfSegment]
    points: list[SrfPoint]


class SrfParseError(Exception):
    pass


PLANE_COUNT_RE = r"PLANE (\d+)"


def read_version(srf_file: TextIO):
    return float(srf_file.readline())


def write_version(srf_file: TextIO):
    srf_file.write("1.0\n")


def read_srf_headers(srf_file: TextIO) -> list[SrfSegment]:
    plane_count_line = srf_file.readline().strip()
    plane_count_match = re.match(PLANE_COUNT_RE, plane_count_line)
    if not plane_count_match:
        raise SrfParseError(f'Expecting PLANE header line, got: "{plane_count_line}"')
    plane_count = int(plane_count_match.group(1))
    segments = []
    for _ in range(plane_count):
        elon = read_float(srf_file)
        elat = read_float(srf_file)
        nstk = read_int(srf_file)
        ndip = read_int(srf_file)
        len = read_float(srf_file)
        wid = read_float(srf_file)
        stk = read_float(srf_file)
        dip = read_float(srf_file)
        dtop = read_float(srf_file)
        shyp = read_float(srf_file)
        dhyp = read_float(srf_file)

        segments.append(
            SrfSegment(
                elon,
                elat,
                nstk,
                ndip,
                len,
                wid,
                stk,
                dip,
                dtop,
                shyp,
                dhyp,
            )
        )
    return segments


POINT_COUNT_RE = r"POINTS (\d+)"


def read_float(srf_file: TextIO, label=None) -> float:
    while (cur := srf_file.read(1)).isspace():
        pass
    float_str = cur
    while not (cur := srf_file.read(1)).isspace():
        float_str += cur
    try:
        return float(float_str)
    except ValueError:
        if label:
            raise SrfParseError(f'Expecting float ({label}), got: "{float_str}"')
        else:
            raise SrfParseError(f'Expecting float, got: "{float_str}"')


def read_int(srf_file: TextIO, label=None) -> int:
    while (cur := srf_file.read(1)).isspace():
        pass
    int_str = cur
    while not (cur := srf_file.read(1)).isspace():
        int_str += cur
    try:
        return int(int_str)
    except ValueError:
        if label:
            raise SrfParseError(f'Expecting int ({label}), got: "{int_str}"')
        else:
            raise SrfParseError(f'Expecting int, got: "{int_str}"')


def read_points_count(srf_file: TextIO) -> int:
    points_count_line = srf_file.readline().strip()
    points_count_match = re.match(POINT_COUNT_RE, points_count_line)
    if not points_count_match:
        raise SrfParseError(f'Expecting POINTS header line, got: "{points_count_line}"')
    return int(points_count_match.group(1))


def read_srf_n_points(point_count: int, srf_file: TextIO) -> list[SrfPoint]:
    points = []
    for _ in range(point_count):
        lon = read_float(srf_file, label="lon")
        lat = read_float(srf_file, label="lat")
        dep = read_float(srf_file, label="dep")
        stk = read_float(srf_file, label="stk")
        dip = read_float(srf_file, label="dip")
        area = read_float(srf_file, label="area")
        tinit = read_float(srf_file, label="tinit")
        dt = read_float(srf_file, label="dt")
        rake = read_float(srf_file, label="rake")
        slip1 = read_float(srf_file, label="slip1")
        nt1 = read_int(srf_file, label="nt1")
        slip2 = read_float(srf_file, label="slip2")
        nt2 = read_int(srf_file, label="nt2")
        slip3 = read_float(srf_file, label="slip3")
        nt3 = read_int(srf_file, label="nt3")
        slipt1 = [read_float(srf_file, label="slipt1") for _ in range(nt1)]
        slipt2 = [read_float(srf_file, label="slipt2") for _ in range(nt2)]
        slipt3 = [read_float(srf_file, label="slipt3") for _ in range(nt3)]
        points.append(
            SrfPoint(
                lon,
                lat,
                dep,
                stk,
                dip,
                area,
                tinit,
                dt,
                rake,
                slip1,
                slip2,
                slip3,
                sr1=slipt1,
                sr2=slipt2,
                sr3=slipt3,
            )
        )
    return points


def read_srf_file(srf_filepath: Path) -> SrfFile:
    with open(srf_filepath, "r") as srf_file:
        srf_file.readline()  # skip version
        header = read_srf_headers(srf_file)
        points = read_srf_points(srf_file)
    return SrfFile(header, points)


def write_srf_header(srf_file: TextIO, header: list[SrfSegment]):
    srf_file.write(f"PLANE {len(header)}\n")
    for segment in header:
        srf_file.write(
            f"{segment.elon:.6f} {segment.elat:.6f} {segment.nstk} {segment.ndip} {segment.len:.4f} {segment.wid:.4f}\n"
        )
        srf_file.write(
            f"{segment.stk:g} {segment.dip:g} {segment.dtop:.4f} {segment.shyp:.4f} {segment.dhyp:.4f}\n"
        )


def write_point_count(srf_file: TextIO, point_count: int):
    srf_file.write(f"POINTS {point_count}\n")


def write_slip(srf_file: TextIO, slips: list[float]):
    for i in range(0, len(slips), 6):
        srf_file.write(
            "  "
            + "  ".join(f"{slips[j]:.5E}" for j in range(i, min(len(slips), i + 6)))
            + "\n"
        )


def write_srf_point(srf_file: TextIO, point: SrfPoint):
    srf_file.write(
        f"{point.lon:.6f} {point.lat:.6f} {point.dep:g} {point.stk:g} {point.dip:g} {point.area:.4E} {point.tinit:.4f} {point.dt:.6E}\n"
    )
    srf_file.write(
        f"{point.rake:g} {point.slip1:.4f} {len(point.sr1)} {point.slip2:.4f} {len(point.sr2)} {point.slip3:.4f} {len(point.sr3)}\n"
    )
    if point.sr1:
        write_slip(srf_file, point.sr1)
    if point.sr2:
        write_slip(srf_file, point.sr2)
    if point.sr3:
        write_slip(srf_file, point.sr3)

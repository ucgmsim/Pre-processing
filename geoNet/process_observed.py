"""Processes V1A or miniseed files from observed events and saves these as text files (in units g)"""
import argparse
from pathlib import Path
from typing import List
import multiprocessing as mp

import numpy as np
import obspy
from scipy import signal
from obspy.clients.fdsn import Client as FDSN_Client

import qcore.timeseries as ts
from geoNet.geoNet_file import GeoNet_File
from geoNet.process import adjust_gf_for_time_delay, SMD

SMD_DEFAULT_KWARGS = {
    "lowcut": 0.05,
    "gpInt": False,
    "highcut": 50.0,
    "fs": (1.0 / 0.005),
    "order": 4,
    "ft": 1.0,
    "output": None,
}

g = 9.810


TEXT_OUTPUT_TYPES = ["accBB", "velLF", "velBB"]


class Record:
    """Represents a observed event record"""

    def __init__(
        self,
        acc_000: np.ndarray,
        acc_090: np.ndarray,
        acc_v: np.ndarray,
        dt: float,
        station_name: str,
        ffp: Path,
    ):
        self.acc_000 = acc_000
        self.acc_090 = acc_090
        self.acc_v = acc_v

        self.rot_angle = None
        self.dt = dt

        self.station_name = station_name
        self.ffp = ffp

    def rotate(self, angle_000: float):
        """
                N
            W<-- -->E
                S

        GeoNet comp_1st and comp_2nd angles axis angles are measured
        from N. We rotate clock wise with angle theta where

            theta = 360 - comp_1st.angle

        [comp_090, comp_180] = Rot(theta) * [comp_1st, comp_2nd]
        Then perform a reflection in the y-axis
        comp_000 = -comp_180
        self.rot_angle:
               theta = 360 - comp_1st.angle
        acc_000:
                positive axis 1 parallel to E
        acc_090:
                positive axis 2 parallel to N

        """
        # For compatibility with orientation as defined by Brendon et al.
        self.rot_angle = -(90.0 - angle_000)
        R = rot_matrix(self.rot_angle)
        acc_090 = R[0, 0] * self.acc_000 + R[0, 1] * self.acc_090
        acc_000 = R[1, 0] * self.acc_000 + R[1, 1] * self.acc_090

        self.acc_090, self.acc_000 = acc_090, acc_000


def pre_process(acc: np.ndarray):
    """De-mean and de-trends the acceleration series"""
    return signal.detrend(acc - acc.mean(), type="linear")


def rot_matrix(theta: float):
    """
    Compute the rotation matrix

    Parameters
    ----------
    theta: float
        Rotation angle in degrees

    Returns
    -------
    array:
        Rotation matrix
    """
    theta = theta * np.pi / 180.0
    return np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])


def process(record: Record, output_format: str, output_type: str = None, **kwargs):
    # Fill any parameter gaps with the default params
    cur_default_kwargs = SMD_DEFAULT_KWARGS.copy()
    cur_default_kwargs["fs"] = 1.0 / record.dt
    cur_default_kwargs["highcut"] = min(1.0 / (2.5 * record.dt), 50)
    kwargs = {**cur_default_kwargs, **kwargs}

    # Perform integration
    smd_comp_000 = SMD(record.acc_000, record.dt, **kwargs)
    smd_comp_090 = SMD(record.acc_090, record.dt, **kwargs)
    smd_comp_v = SMD(record.acc_v, record.dt, **kwargs)

    if output_format == "text":
        assert output_type in TEXT_OUTPUT_TYPES
        output_dir = record.ffp.parent / output_type
        output_dir.mkdir(exist_ok=True)
        ts.seis2txt(
            smd_comp_000.__getattribute__(output_type),
            smd_comp_000.dt,
            str(output_dir) + "/",
            record.station_name,
            "000",
        )
        ts.seis2txt(
            smd_comp_090.__getattribute__(output_type),
            smd_comp_090.dt,
            str(output_dir) + "/",
            record.station_name,
            "090",
        )
        ts.seis2txt(
            smd_comp_v.__getattribute__(output_type),
            smd_comp_v.dt,
            str(output_dir) + "/",
            record.station_name,
            "ver",
        )
    else:
        raise NotImplementedError()


def get_comp_trace(traces: List[obspy.Trace], component: str):
    """Retrieves the trace for the specified component"""
    return [
        cur_trace.data
        for cur_trace in traces
        if cur_trace.stats.channel[-1] == component
    ][0]


def process_v1a_file(ffp: Path, output_format: str, output_type, ix: int, n_files: int):
    """Processes a single V1A file and saves as text file, to be used with pool.starmap"""
    print(f"Processing V1A file {ffp.name}, {ix + 1}/{n_files}")
    try:
        gf = GeoNet_File(ffp.name, base_dir=str(ffp.parent))

        if gf.comp_1st.acc.size < 5.0 / gf.comp_1st.delta_t:
            print(f"{ffp.name} has less than 5 secs of data, skipping!")
            return None

        # When appended zeroes at the beginning of the record are removed, the
        # record might then be empty, skip processing in such a case
        agf = adjust_gf_for_time_delay(gf)
        if agf.comp_1st.acc.size <= 10:
            print(f"No elements in {ffp.name}, skipping!")
            return None

        # Pre-processing (de-mean & detrend)
        acc_000, acc_090 = (
            pre_process(gf.comp_1st.acc),
            pre_process(gf.comp_2nd.acc),
        )
        acc_v = pre_process(gf.comp_up.acc)

        # Create a record and rotate
        cur_record = Record(
            acc_000,
            acc_090,
            acc_v,
            gf.comp_1st.delta_t,
            ffp.name.split(".")[0].split("_")[2],
            ffp,
        )
        cur_record.rotate(gf.comp_1st.angle)

        # Process
        process(cur_record, output_format, output_type)
        return True

    except Exception as ex:
        print(f"Failed to process V1A file {ffp} with exception:\n {ex}")
        return False


def process_v1a_files(
    ffps: List[Path], output_format: str, output_type: str = None, n_procs: int = 1
):
    """Processes the specified V1A files"""
    n_files = len(ffps)
    with mp.Pool(n_procs) as p:
        success_mask = p.starmap(
            process_v1a_file,
            [
                (cur_ffp, output_format, output_type, ix, n_files)
                for ix, cur_ffp in enumerate(ffps)
            ],
        )

    return success_mask


def process_miniseed_file(
    ffp: Path,
    inventory: obspy.Inventory,
    output_format: str,
    output_type: str,
    ix: int,
    n_files: int,
):
    """Processes a single miniseed file and saves as text file, to be used with pool.starmap"""
    print(f"Processing miniseed file {ffp.name}, {ix + 1}/{n_files}")
    try:
        st = obspy.read(str(ffp))
        if len(st) == 3:
            # Remove sensitivity
            st = st.remove_sensitivity(inventory=inventory)

            # Rotate if required
            st.rotate("->ZNE", inventory=inventory)

            # Convert to ACC if it isn't already
            assert st.traces[0].stats.channel[-2] in ["N", "H"]
            if st.traces[0].stats.channel[-2] == "H":
                st.differentiate()

            # Pre-processing (de-mean & detrend) & convert to units g
            acc_000 = pre_process(get_comp_trace(st.traces, "N")) / g
            acc_090 = pre_process(get_comp_trace(st.traces, "E")) / g
            acc_v = pre_process(get_comp_trace(st.traces, "Z")) / g

            # Process
            cur_record = Record(
                acc_000,
                acc_090,
                acc_v,
                st.traces[0].stats.delta,
                st.traces[0].stats.station,
                ffp,
            )
            process(cur_record, output_format, output_type)
            return True
        else:
            print(f"File {ffp.name} does not have 3 components, skipping!")
            return None
    except Exception as ex:
        print(f"Failed to process miniseed file {ffp} with exception:\n {ex}")
        return False


def process_miniseed_files(
    ffps: List[Path], output_format: str, output_type: str = None, n_procs: int = 1
):
    """Processes the specified miniseed files"""
    print("Loading the station inventory (this may take a few seconds)")
    client = FDSN_Client("GEONET")
    inventory = client.get_stations(station="*", level="response")

    n_files = len(ffps)
    with mp.Pool(n_procs) as p:
        success_mask = p.starmap(
            process_miniseed_file,
            [
                (cur_ffp, inventory, output_format, output_type, ix, n_files)
                for ix, cur_ffp in enumerate(ffps)
            ],
        )

    return success_mask


def main(
    base_data_dir: Path,
    input_format: str,
    output_format: str,
    output_type: str,
    n_procs: int,
):
    # Find all files to process
    print("Finding all relevant files to process")
    ffps = list(base_data_dir.rglob(f"*.{input_format}"))
    print(f"Found {len(ffps)} files to process")

    if input_format.lower() == "v1a":
        success_mask = process_v1a_files(
            ffps, output_format, output_type=output_type, n_procs=n_procs
        )
    else:
        success_mask = process_miniseed_files(
            ffps, output_format, output_type=output_type, n_procs=n_procs
        )

    unprocessed_files = [
        f"{ffps[ix]}\n" for ix, result in enumerate(success_mask) if result is None
    ]
    if len(unprocessed_files) > 0:
        print(
            f"The following were not processed as they did not meet the data requirements."
            f"Check the log for the reason.\n{''.join(unprocessed_files)}"
        )

    failed_files = [
        f"{ffps[ix]}\n" for ix, result in enumerate(success_mask) if result is False
    ]
    if len(failed_files) > 0:
        print(
            f"The following files failed to process, "
            f"check the log for exceptions.\n{''.join(failed_files)}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Converts V1A or miniseed files to processed binary or text files."
    )

    parser.add_argument(
        "base_data_dir",
        type=str,
        help="The base directory for the observed GM records. "
        "All V1A or miniseed files below this directory will be processed.",
    )
    parser.add_argument(
        "input_format",
        type=str,
        help="The data format of the observed record files",
        choices=["V1A", "mseed"],
    )
    parser.add_argument(
        "output_format",
        type=str,
        help="The format of the output file, either 3 text files "
        "(one for each component) or a single binary file.",
        choices=["text", "binary"],
    )
    parser.add_argument(
        "--output_type",
        type=str,
        help="Output type to save the data in, only relevant when saving in text format",
        default="accBB",
        choices=TEXT_OUTPUT_TYPES,
    )
    parser.add_argument(
        "--n_procs", type=int, help="Number of processes to use", default=1
    )

    args = parser.parse_args()
    main(
        Path(args.base_data_dir),
        args.input_format,
        args.output_format,
        args.output_type,
        args.n_procs,
    )

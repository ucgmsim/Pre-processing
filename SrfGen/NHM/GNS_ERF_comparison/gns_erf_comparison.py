import os

import numpy as np

from qcore.nhm import load_nhm, NHMFault


class GNSFault:
    def __init__(
        self,
        fault_name,
        slip_type,
        n_line_sections,
        dip,
        dip_dir,
        depth_to_base,
        depth_to_top,
        end_1,
        end_2,
        mag,
        rp,
        line_segments,
        trace,
    ):
        self.fault_name = fault_name
        self.slip_type = slip_type
        self.n_line_sections = n_line_sections
        self.dip = dip
        self.dip_dir = dip_dir
        self.depth_to_base = depth_to_base
        self.depth_to_top = depth_to_top
        self.end_1 = end_1
        self.end_2 = end_2
        self.mag = mag
        self.rp = rp
        self.line_segments = line_segments
        self.trace = trace

    @staticmethod
    def read_fault(lines, line_ix):
        def read_loc_line(cur_line):
            cur_line = [float(entry) for entry in cur_line.split()]
            loc_1 = (cur_line[2] + (cur_line[3] / 60), -1 * (cur_line[0] + (cur_line[1] / 60)))
            loc_2 = (cur_line[6] + (cur_line[7] / 60), -1 * (cur_line[4] + (cur_line[5] / 60)))

            return loc_1, loc_2

        # Greendale     ss
        fault_name, slip_type = lines[line_ix].split()

        # 6D
        n_line_sections = int(lines[line_ix + 1].strip()[:-1])

        # 85.0 180.0 10.0  0.0
        dip, direction, depth_to_base, depth_to_top = [
            float(entry) for entry in lines[line_ix + 2].split()
        ]

        # 43 34.6 172  1.8  43 34.5 172 22.5         7.1                  25000
        cur_line = [float(entry) for entry in lines[line_ix + 3].split()]
        end_1 = (cur_line[2] + (cur_line[3] / 60), cur_line[0] + (cur_line[1] / 60))
        end_2 = (cur_line[6] + (cur_line[7] / 60), cur_line[4] + (cur_line[5] / 60))
        mag, rp = cur_line[8], cur_line[9]

        line_segments = []
        for ix in range(n_line_sections):
            line_segments.append(read_loc_line(lines[line_ix + 4 + ix]))

        trace_points = []
        for ix in range(len(line_segments)):
            gns_loc_1, gns_loc_2 = line_segments[ix]
            trace_points.append(gns_loc_1)
            if ix + 1 == len(line_segments):
                trace_points.append(gns_loc_2)

        trace_points = np.array(trace_points)

        return (
            GNSFault(
                fault_name,
                slip_type,
                n_line_sections,
                dip,
                direction,
                depth_to_base,
                depth_to_top,
                end_1,
                end_2,
                mag,
                rp,
                line_segments,
                trace_points
            ),
            line_ix + 4 + n_line_sections + 1,
        )


def compare(gns_fault: GNSFault, nhm_fault: NHMFault):
    result_str, diff = f"\n{nhm_fault.name} has different:\n", False
    if nhm_fault.dip != gns_fault.dip:
        result_str += f"\t - Dip: GNS - {gns_fault.dip}, NHM - {nhm_fault.dip}\n"
        diff = True
    if not np.isclose(nhm_fault.dip_dir, gns_fault.dip_dir, rtol=0.01):
        result_str += (
            f"\t - Dip direction: GNS - {gns_fault.dip_dir}, NHM - {nhm_fault.dip_dir}\n"
        )
        diff = True
    # if np.round(nhm_fault.mw, 1) != gns_fault.mag:
    if not np.isclose(nhm_fault.mw, gns_fault.mag, 0.01):
        result_str += f"\t - Magnitude: GNS - {gns_fault.mag}, NHM - {nhm_fault.mw}\n"
        diff = True
    # if (
    #     (np.round(nhm_fault.recur_int_median, -1) != gns_fault.rp)
    #     and (np.round(nhm_fault.recur_int_median, -2) != gns_fault.rp)
    #     and (np.round(nhm_fault.recur_int_median, -3) != gns_fault.rp)
    #     and (np.round(nhm_fault.recur_int_median, -4) != gns_fault.rp)
    # ):
    if not np.isclose(nhm_fault.recur_int_median, gns_fault.rp, rtol=0.05):
        result_str += f"\t - Return period: GNS - {gns_fault.rp}, NHM - {nhm_fault.recur_int_median}\n"
        diff = True
    if nhm_fault.dtop != gns_fault.depth_to_top:
        result_str += f"\t - Dtop: GNS - {gns_fault.depth_to_top}, NHM - {nhm_fault.dtop}\n"
        diff = True
    if nhm_fault.dbottom != gns_fault.depth_to_base:
        result_str += f"\t - Dbottom: GNS - {gns_fault.depth_to_base}, NHM - {nhm_fault.dbottom}\n"
        diff = True

    # NHM contains points along the trace, whereas GNS is individual line segments
    if gns_fault.trace.shape[0] != nhm_fault.trace.shape[0] or not np.all(np.isclose(nhm_fault.trace, gns_fault.trace)):
        result_str += f"\t - Fault traces:\nGNS \n{gns_fault.trace}\nNHM \n{nhm_fault.trace}\n"
        diff = True

    # # Compare line segments
    # for ix in range(nhm_fault.trace.shape[0]):
    #     nhm_loc = nhm_fault.trace[ix]
    #     if (ix + 1) < nhm_fault.trace.shape[0]:
    #         gns_loc_1, gns_loc_2 = gns_fault.line_segments[ix]
    #         gns_loc = np.asarray(gns_loc_1)
    #     else:
    #         gns_loc_1, gns_loc_2 = gns_fault.line_segments[ix - 1]
    #         gns_loc = np.asarray(gns_loc_2)
    #
    #     if not np.all(np.isclose(nhm_loc, gns_loc)):
    #         diff = True

    if diff:
        return result_str
    return None


if __name__ == "__main__":
    standard_erf_ffp = (
        os.path.joins(os.path.dirname(__file__), "../NZ_FLTmodel_2010_v18p6.txt")
    )
    gns_erf_ffp = "./F501111U.DAT"

    out_file = ""
    if len(out_file) == 0 or os.path.isfile(out_file):
        print("Specify a valid output file, quitting.")
        exit()

    nhm_faults = load_nhm(standard_erf_ffp)

    with open(gns_erf_ffp) as f:
        lines = f.readlines()

    # Ignore the first 3 lines
    lines = lines[3:]

    line_ix, gns_faults = 0, []
    while line_ix < len(lines):
        print(f"Processing line {line_ix}/{len(lines)}")
        cur_fault, line_ix = GNSFault.read_fault(lines, line_ix)
        gns_faults.append(cur_fault)

    gns_fault_names = set(
        [cur_fault.fault_name.replace("&", "-") for cur_fault in gns_faults]
    )
    nhm_fault_names = set(list(nhm_faults.keys()))

    shared_names = gns_fault_names.intersection(nhm_faults)
    nhm_missing = gns_fault_names.difference(shared_names)
    gns_missing = nhm_fault_names.difference(shared_names)

    print(f"NHM missing {nhm_missing}")
    print(f"GNS missing {gns_missing}")

    result_out = []
    for cur_gns_fault in gns_faults:
        if cur_gns_fault.fault_name in shared_names:
            cur_result = compare(cur_gns_fault, nhm_faults[cur_gns_fault.fault_name])
            if cur_result:
                result_out.append(cur_result)
                print(cur_result)

    if len(result_out) > 0:
        with open(out_file, 'w') as f:
            f.write(f"Faults not in NHM {nhm_missing}\n\n")
            # f.write(f"Faults not in GNS {gns_missing}\n")
            f.write("Comparison of faults in both:\n")
            f.writelines(result_out)

    exit()

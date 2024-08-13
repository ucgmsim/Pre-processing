"""
readGP and writeGP have been moved to utils. 
read_geoNet_list moved to geoNet_file to make it independent. 
utilities should not be used and will be removed in the future
"""
import numpy as np


def readGP(loc, fname):
    """
    Convinience function for reading files in the Graves and Pitarka format
    """
    with open("/".join([loc, fname]), "r") as f:
        lines = f.readlines()

    data = []

    for line in lines[2:]:
        data.append([float(val) for val in line.split()])

    data = np.concatenate(data)

    return data


def read_geoNet_list(lines, line_width=80, width=8):
    """
    Convinience function for parsing lines in GeoNet format
    """
    data = []
    slices = np.arange(0, line_width, width)
    for line in lines[:-1]:
        for i in slices:
            if line[i : i + width] == "999999.9":
                return np.asarray(data, dtype=float)
            else:
                data.append(float(line[i : i + width]))

    last_line = lines[-1].rstrip()
    for i in xrange(0, len(last_line), width):
        if last_line[i : i + width] == "999999.9":
            return np.asarray(data, dtype=float)
        else:
            data.append(float(last_line[i : i + width]))

    return np.asarray(data, dtype=float)


def writeGP(loc, fname, data, header, ncol=6):
    """
    Convinience function for writing files in the Graves and Pitarka format
    """
    size = len(data)
    nrow = size / ncol
    size_last_row = size % ncol

    lines = ""
    for line in np.reshape(xrange(nrow * ncol), (nrow, ncol)):
        for val in line:
            lines += "{:.6e}".format(data[val]) + 3 * " "
        lines = lines.rstrip(3 * " ") + "\n"

    if size_last_row:
        for i in xrange(1, size_last_row + 1):
            lines += "{:.6e}".format(data[-i]) + 3 * " "
    lines = lines.rstrip(3 * " ")

    with open("/".join([loc, fname]), "w") as f:
        f.writelines(header)
        f.writelines(lines)
    return

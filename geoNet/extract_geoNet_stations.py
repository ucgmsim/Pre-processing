"""
See networks.csv https://github.com/GeoNet/delta/blob/master/network/networks.csv
for what the codes mean. The relevant codes for QuakeCoRE are (I am guessing here)
    NZ (included per Brendon's advice
    SM National strong motion network
    SC Canterbury regional strong motion network
    SB is included in http://info.geonet.org.nz/display/equip/Network+Location+Queries
       operational strong motions stations XML file.
    
    SX private sites strong motion.

stations.csv can be found at
    https://github.com/GeoNet/delta/tree/master/network
"""


import os
import csv
import numpy as np


def In_southIsland(lon, lat):
    """
    define a rectangular box given a centre, coordinates of length and width of the box
    (xp, yp) define the new coordinate system
    """
    centre_lat = -43.894316
    centre_lon = 170.650197

    L_lat = -40.877260
    L_lon = 174.528072

    W_lat = -42.733521
    W_lon = 169.254083

    theta = -45.0

    x0 = centre_lon
    y0 = centre_lat
    x = lon
    y = lat

    xp = np.cos(theta * np.pi / 180.0) * (x - x0) - np.sin(theta * np.pi / 180.0) * (
        y - y0
    )
    yp = np.sin(theta * np.pi / 180.0) * (x - x0) + np.cos(theta * np.pi / 180.0) * (
        y - y0
    )

    xp_L = np.cos(theta * np.pi / 180.0) * (L_lon - x0) - np.sin(
        theta * np.pi / 180.0
    ) * (L_lat - y0)
    yp_L = np.sin(theta * np.pi / 180.0) * (L_lon - x0) + np.cos(
        theta * np.pi / 180.0
    ) * (L_lat - y0)
    L = np.sqrt(xp_L**2 + yp_L**2)

    xp_W = np.cos(theta * np.pi / 180.0) * (W_lon - x0) - np.sin(
        theta * np.pi / 180.0
    ) * (W_lat - y0)
    yp_W = np.sin(theta * np.pi / 180.0) * (W_lon - x0) + np.cos(
        theta * np.pi / 180.0
    ) * (W_lat - y0)
    W = np.sqrt(xp_W**2 + yp_W**2)

    if np.abs(xp) <= L and np.abs(yp) <= W:
        within_southIsland = True
    else:
        within_southIsland = False

    return within_southIsland


def In_northIsland(lon, lat):
    """
    define a rectangular box given a centre, coordinates of length and width of the box
    (xp, yp) define the new coordinate system
    """
    centre_lat = -37.651034
    centre_lon = 175.435162

    L_lat = -34.104810
    L_lon = 172.823488

    W_lat = -38.888624
    W_lon = 172.315141

    theta = -329.0

    x0 = centre_lon
    y0 = centre_lat
    x = lon
    y = lat

    xp = np.cos(theta * np.pi / 180.0) * (x - x0) - np.sin(theta * np.pi / 180.0) * (
        y - y0
    )
    yp = np.sin(theta * np.pi / 180.0) * (x - x0) + np.cos(theta * np.pi / 180.0) * (
        y - y0
    )

    xp_L = np.cos(theta * np.pi / 180.0) * (L_lon - x0) - np.sin(
        theta * np.pi / 180.0
    ) * (L_lat - y0)
    yp_L = np.sin(theta * np.pi / 180.0) * (L_lon - x0) + np.cos(
        theta * np.pi / 180.0
    ) * (L_lat - y0)
    L = np.sqrt(xp_L**2 + yp_L**2)

    xp_W = np.cos(theta * np.pi / 180.0) * (W_lon - x0) - np.sin(
        theta * np.pi / 180.0
    ) * (W_lat - y0)
    yp_W = np.sin(theta * np.pi / 180.0) * (W_lon - x0) + np.cos(
        theta * np.pi / 180.0
    ) * (W_lat - y0)
    W = np.sqrt(xp_W**2 + yp_W**2)

    if np.abs(xp) <= L and np.abs(yp) <= W:
        within_southIsland = True
    else:
        within_southIsland = False

    return within_southIsland


fname = "stations.csv"
f = open("/".join([os.getcwd(), fname]), "r")
fcsv = csv.DictReader(f, delimiter=",")
with open("all_geoNet_stats.ll", "w") as fstats:
    for line in fcsv:
        # Skip stations that have closed
        if line["End Date"] != "9999-01-01T00:00:00Z":
            continue

        Network = line["Network"]
        if Network in ["NZ", "SM", "SC", "SB", "SX"]:
            pass  # don't do nothing
        else:
            continue  # skip loop

        lon = float(line["Longitude"])
        lat = float(line["Latitude"])
        # only save if within Main Land NZ
        if In_southIsland(lon, lat) or In_northIsland(lon, lat):
            pass
        else:
            continue

        fstats.write("%10.4f  %10.4f  %10s\n" % (lon, lat, line["Station"]))

f.close()

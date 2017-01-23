import numpy as np
import os

google_earth_start_kml="""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
    <name>station_locations.kml</name>
    <Style id="sh_ylw-pushpin00">
        <IconStyle>
            <scale>1.3</scale>
            <Icon>
                <href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
            </Icon>
            <hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
        </IconStyle>
        <LineStyle>
            <color>ff0055ff</color>
            <width>3.1</width>
        </LineStyle>
    </Style>
    <Style id="sn_ylw-pushpin21">
        <IconStyle>
            <scale>1.1</scale>
            <Icon>
                <href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
            </Icon>
            <hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
        </IconStyle>
        <LineStyle>
            <color>ff0055ff</color>
            <width>3.1</width>
        </LineStyle>
    </Style>
    <StyleMap id="msn_ylw-pushpin1">
        <Pair>
            <key>normal</key>
            <styleUrl>#sn_ylw-pushpin21</styleUrl>
        </Pair>
        <Pair>
            <key>highlight</key>
            <styleUrl>#sh_ylw-pushpin00</styleUrl>
        </Pair>
    </StyleMap>
"""

google_earth_end_kml="""
</Document>
</kml>
"""

#def read_statsll(fname=None, loc=os.getcwd()):
#    stats_dict = {}
#    with open("/".join([loc, fname]), "r") as f:
#        for line in f:
#            lon, lat, stat_code = line.split()
#            stat_code = stat_code.strip()
#            stats_dict.update({stat_code: [float(lon), float(lat)]})
#
#    return stats_dict
#
#statsll_fname = "20150105_175001_eventStats_2016-12-12.ll" 
#stats_dict = read_statsll(fname=statsll_fname,
#                         loc=".")
#
#with open("20150105_175001_eventStats_2016-12-12.kml", 'w') as f:
#    f.write(google_earth_start_kml)
#    for stat_code in stats_dict.keys():
#        f.write("\n")
#        f.write("<Placemark>")
#        f.write("<name> %s </name>" %stat_code)
#        f.write("<styleUrl>#msn_ylw-pushpin1</styleUrl>")
#        f.write("<Point>")
#        f.write("<gx:drawOrder>1</gx:drawOrder>")
#        lon, lat = stats_dict[stat_code]
#        f.write("<coordinates>%s,%s,0</coordinates>" %(str(lon), str(lat)))
#        f.write("</Point>")
#        f.write("</Placemark>")
#        f.write("\n")
#    f.write(google_earth_end_kml)

def write_stats_kml(loc, fname, stats_dict):
    """
    stats_dict: dict in the form given by read_statsll
    """
    if not fname.endswith(".kml"):
        fname+=".kml"

    with open("/".join([loc, fname]), 'w') as f:
        f.write(google_earth_start_kml)
        for stat_code in stats_dict.keys():
            f.write("\n")
            f.write("<Placemark>")
            f.write("<name> %s </name>" %stat_code)
            f.write("<styleUrl>#msn_ylw-pushpin1</styleUrl>")
            f.write("<Point>")
            f.write("<gx:drawOrder>1</gx:drawOrder>")
            lon, lat = stats_dict[stat_code]
            f.write("<coordinates>%s,%s,0</coordinates>" %(str(lon), str(lat)))
            f.write("</Point>")
            f.write("</Placemark>")
            f.write("\n")
        f.write(google_earth_end_kml)
    return

def gen_stats_kml(loc,fname):
    """
    fname: station_file_name.ll
    saves: station_file_name.kml
    """
    return

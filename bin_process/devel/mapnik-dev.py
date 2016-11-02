#!/usr/bin/env python2

import mapnik as mn

# Projection Definitions - spatialreference.org
GOOGLE = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs'
NZGD2000 = '+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
LONGLAT = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'

# set projection of displayed map
MAP_PROJECTION = GOOGLE

def proj_transform(proj_src, proj_dest, location):
    # these have to be variables for C++ backend to not segfault
    a = mn.Projection(proj_src)
    b = mn.Projection(proj_dest)
    return mn.ProjTransform(a, b).forward(location)

# map of length, width pixels
m = mn.Map(3400, 4000)
# set projection
m.srs = MAP_PROJECTION
# fill
m.background = mn.Color('lightblue')

s = mn.Style()
r = mn.Rule()
pen = mn.LineSymbolizer()
pen.stroke = mn.Color('rgb(100%,100%,100%)')
pen.stroke_width = 1.0
r.symbols.append(pen)
s.rules.append(r)
m.append_style('white line',s)

s = mn.Style()
r = mn.Rule()
pen = mn.LineSymbolizer()
pen.stroke = mn.Color('yellow')
pen.stroke_width = 2.0
r.symbols.append(pen)
s.rules.append(r)
m.append_style('shwy line',s)

s = mn.Style()
r = mn.Rule()
filler = mn.PolygonSymbolizer()
filler.fill = mn.Color('lightblue')
r.symbols.append(filler)
s.rules.append(r)
m.append_style('water',s)

s = mn.Style()
r = mn.Rule()
filler = mn.PolygonSymbolizer()
filler.fill = mn.Color('rgb(120, 150, 120)')
r.symbols.append(filler)
pen = mn.LineSymbolizer()
pen.stroke = mn.Color('rgb(0%,0%,0%)')
pen.stroke_width = 1.0
r.symbols.append(pen)
s.rules.append(r)
m.append_style('land fill',s)

highways = mn.Layer('highways', NZGD2000)
highways.datasource = mn.Shapefile(file='/home/vap30/shared/geo/shwy/shwy.shp')
highways.styles.append('shwy line')

roads = mn.Layer('roads', LONGLAT)
roads.datasource = mn.Shapefile(file='/home/vap30/shared/geo/lds-nz-road-centre-line/nz-road-centre-line-electoral.shp')
roads.styles.append('white line')

land = mn.Layer('land', LONGLAT)
land.datasource = mn.Shapefile(file='/home/vap30/shared/geo/lds-nz-coastlines-and-islands/nz-coastlines-and-islands-polygons-topo-1500k.shp')
land.styles.append('land fill')

lakes = mn.Layer('lakes', LONGLAT)
lakes.datasource = mn.Shapefile(file='/home/vap30/shared/geo/lds-nz-lake-polygons/nz-lake-polygons-topo-1500k.shp')
lakes.styles.append('water')

m.layers.append(land)
m.layers.append(lakes)
m.layers.append(roads)
m.layers.append(highways)

# crop to region
ll_region = mn.Envelope(166, -47.5, 174.5, -40)
bbox = proj_transform(LONGLAT, MAP_PROJECTION, ll_region)
m.zoom_to_box(bbox)
# render
mn.render_to_file(m, 'world.svg')

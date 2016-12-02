"""
Library which simplifies the use of GMT.
Functions should hide obvious parameters.
"""

import os
from shutil import copyfile, move
from subprocess import call, PIPE, Popen
from sys import byteorder
from time import time

import numpy as np

# only needed if plotting fault planes direct from SRF
try:
    from shared_srf import *
except ImportError:
    print('shared_srf.py not found. will not be able to plot faults from SRF.')

# if gmt available in $PATH, GMT should be 'gmt'
# to use a custom location, set full path to 'gmt' binary below
GMT = 'gmt'

# retrieve version of GMT
gmtp = Popen([GMT, '--version'], stdout = PIPE)
GMT_VERSION = gmtp.communicate()[0].rstrip()
GMT_MAJOR, GMT_MINOR = map(int, GMT_VERSION.split('.')[:2])
if GMT_MAJOR != 5:
    print('This library is only for GMT version 5. You have %s.' \
            % (GMT_VERSION))
# ps2raster becomes psconvert in GMT 5.2
elif GMT_MINOR < 2:
    psconvert = 'ps2raster'
else:
    print('WARNING: LIBRARY NOT TESTED WITH GMT 5.2 CHANGES')
    print('WILL MOST DEFINITLY CAUSE ISSUES, OR USE 5.1')
    psconvert = 'psconvert'

###
### COMMON RESOURCES
###
# definition of locations which can be mapped
# longitude, latitude,
# point position [Left Centre Right, Top Middle Bottom]
sites = { \
    'Akaroa':(172.9683333, -43.80361111, 'RB'), \
    'Blenheim':(173.9569444, -41.5138888, 'LM'), \
    'Christchurch':(172.6347222, -43.5313888, 'LM'), \
    'Darfield':(172.1116667, -43.48972222, 'CB'), \
    'Dunedin':(170.3794444, -45.8644444, 'LM'), \
    'Greymouth':(171.2063889, -42.4502777, 'RM'), \
    'Haast':(169.0405556, -43.8808333, 'LM'), \
    'Kaikoura':(173.6802778, -42.4038888, 'LM'), \
    'Lyttleton':(172.7194444, -43.60305556, 'LM'), \
    'Masterton':(175.658333, -40.952778, 'LM'), \
    'Napier':(176.916667, -39.483333, 'LM'), \
    'New Plymouth':(174.083333, -39.066667, 'RM'), \
    'Nelson':(173.2838889, -41.2761111, 'CB'), \
    'Oxford':(172.1938889, -43.29555556, 'LB'), \
    'Palmerston North':(175.611667, -40.355000, 'RM'), \
    'Queenstown':(168.6680556, -45.0300000, 'LM'), \
    'Rakaia':(172.0230556, -43.75611111, 'RT'), \
    'Rolleston':(172.3791667, -43.59083333, 'RB'), \
    'Rotorua':(176.251389, -38.137778, 'LM'), \
    'Taupo':(176.069400, -38.6875, 'LM'), \
    'Tekapo':(170.4794444, -44.0069444, 'LM'), \
    'Timaru':(171.2430556, -44.3958333, 'LM'), \
    'Wellington':(174.777222, -41.288889, 'RM'), \
    'Westport':(171.5997222, -41.7575000, 'RM')}

###
### ACCESSORY FUNCTIONS
###
def is_native_xyv(xyv_file, x_min, x_max, y_min, y_max, v_min = None):
    """
    Detects whether an input file is native or if it needs bytes swapped.
    It makes sure values are sane, if not. Non-native is assumed.
    xyv_file: file containing 3 columns
    x_min: minimum x value (first column)
    x_max: maximum x value
    y_min: minimum y value (second column)
    y_max: maximum y value
    v_min: minimum value in third column (None for skip)
    """
    # form array of xyv data (3 columns of 4 byte floats)
    bin_data = np.fromfile(xyv_file, dtype = '3f4')

    # check the first few rows
    for i in xrange(min(10, len(bin_data))):
        if x_min <= bin_data[i, 0] <= x_max and \
                y_min <= bin_data[i, 1] <= y_max and \
                (v_min == None or v_min <= bin_data[i, 2]):
            continue
        else:
            # invalid values, not native
            return False
    # no invalid values found, assuming native endian
    return True

def swap_bytes(xyv_file, native_version, bytes_per_var = 4):
    """
    Simple and fast way to swap bytes in a file.
    xyv_file: input file
    native_version: where to store the result
    bytes_per_var: how long each value is
    """
    if byteorder == 'little':
        data = np.fromfile(xyv_file, dtype = '>f%d' % (bytes_per_var))
    else:
        data = np.fromfile(xyv_file, dtype = '<f%d' % (bytes_per_var))

    data.tofile(native_version)

def abs_max(x_file, y_file, z_file, out_file, native = True):
    """
    Creates a file containing the absolute max value of 3 components.
    Each file is assumed to contain 3 columns of 4 byte values.
    x_file: 1st input file (named x here)
    y_file: 2nd input file (named y here)
    z_file: 3rd input file (named z here)
    out_file: where to store the result
    native: files are in native endian if True
    """
    # allow all-in-one with byteswap capability
    if native:
        fmt = '3f4'
    elif byteorder == 'little':
        fmt = '3>f4'
    else:
        fmt = '3<f4'

    result = np.fromfile(x_file, dtype = fmt)
    y = np.fromfile(y_file, dtype = fmt)[:, 2]
    z = np.fromfile(z_file, dtype = fmt)[:, 2]

    result[:, 2] = np.sqrt(result[:, 2] ** 2 + y ** 2 + z ** 2)
    result.astype('f4').tofile(out_file)

# TODO: function should be able to modify result CPT such that:
#       background colour is extended just like foreground (bidirectional)
def makecpt(source, output, low, high, inc = None, invert = False):
    """
    Creates a colour palette file.
    source: inbuilt scale or template file
    output: filepath to store file
    low: minimum range
    high: maximum range
    inc: discrete increment
    invert: whether to swap colour order
    """
    # work out GMT colour range parameter
    crange = '%s/%s' % (low, high)
    if inc != None:
        crange = '%s/%s' % (crange, inc)

    cmd = [GMT, 'makecpt', '-T%s' % (crange), '-C%s' % (source)]
    if invert:
        cmd.append('-I')
    with open(output, 'w') as cptf:
        Popen(cmd, stdout = cptf).wait()

def grd_mask(xy_file, out_file, region = None, dx = '1k', dy = '1k'):
    """
    Creates a mask file from a path. Inside = 1, outside = NaN.
    xy_file: file containing a path
    out_file: name of output GMT grd file
    region: tuple region of grd file (must be set if gmt.history doesn't exist)
    dx: x grid spacing size
    dy: y grid spacing size
    """
    wd = os.path.dirname(out_file)
    out_file = os.path.basename(out_file)
    if wd == '':
        wd = '.'
    cmd = ([GMT, 'grdmask', xy_file, '-G%s' % (out_file), '-NNaN/1/1', \
            '-I%s/%s' % (dx, dy)])
    if region == None:
        cmd.append('-R')
    else:
        cmd.append('-R%s/%s/%s/%s' % region)
    Popen(cmd, cwd = wd).wait()

def gmt_defaults(wd = '.', font_annot_primary = 16, map_tick_length_primary = '0.05i', \
        font_label = 16, ps_page_orientation = 'portrait', map_frame_pen = '1p', \
        format_geo_map = 'D', map_frame_type = 'plain', format_float_out = '%lg', \
        proj_length_unit = 'i', ps_media = 'A2', extra = []):
    """
    Sets default values for GMT.
    GMT stores these values in the file 'gmt.conf'
    wd: which directory to set for
    extra: list of params eg: ['FONT_ANNOT_SECONDARY', '12', 'KEY', '=', 'VALUE']
    """
    cmd = ['gmtset', \
            'FONT_ANNOT_PRIMARY', '%s' % (font_annot_primary), \
            'MAP_TICK_LENGTH_PRIMARY', '%s' % (map_tick_length_primary), \
            'FONT_LABEL', '%s' % (font_label), \
            'PS_PAGE_ORIENTATION', ps_page_orientation, \
            'MAP_FRAME_PEN', '%s' % (map_frame_pen), \
            'FORMAT_GEO_MAP', format_geo_map, \
            'MAP_FRAME_TYPE', map_frame_type, \
            'FORMAT_FLOAT_OUT', format_float_out, \
            'PROJ_LENGTH_UNIT', proj_length_unit, \
            'PS_MEDIA', '=', ps_media]
    # protect users from entering non-string values
    cmd.extend(map(str, extra))
    Popen(cmd, cwd = wd).wait()

###
### MAIN PLOTTING CLASS
###
class GMTPlot:

    def __init__(self, pspath, append = False):
        self.pspath = pspath
        if append:
            self.psf = open(pspath, 'a')
            self.new = False
        else:
            self.psf = open(pspath, 'w')
            self.new = True
        # figure out where to run GMT from
        self.wd = os.path.abspath(os.path.dirname(pspath))
        if self.wd == '':
            self.wd = os.path.abspath('.')
        # if previous/custom gmt.defaults and gmt.history was wanted
        # it should have already been copied into this directory
        if not os.path.exists(os.path.join(self.wd, 'gmt.defaults')):
            gmt_defaults(wd = self.wd)
        # place to reject unwanted warnings
        self.sink = open('/dev/null', 'a')

    def background(self, length, height, \
            left_margin = 0, bottom_margin = 0, colour = 'white'):
        """
        Draws background on GMT plot.
        This should be the first action.
        length: how wide the background should be
        height: how high the background should be
        left_margin: to move origin left afterwards
        bottom_margin: to move origin up afterwards
        colour: the colour of the background
        """
        # draw background and place origin up, right as wanted
        cmd = [GMT, 'psxy', '-K', '-G%s' % (colour), \
                '-JX%s/%s' % (length, height), '-R0/%s/0/%s' % (length, height), \
                '-Xa-%s' % (left_margin), '-Ya-%s' % (bottom_margin)]
        # one of the functions that can be run on a blank file
        # as such, '-O' flag needs to be taken care of
        if self.new:
            self.new = False
        else:
            cmd.append('-O')
        proc = Popen(cmd, stdin = PIPE, stdout = self.psf, cwd = self.wd)
        proc.communicate('%s 0\n%s %s\n0 %s\n0 0' \
                % (length, length, height, height))
        proc.wait()

    def spacial(self, proj, region, \
            lon0 = None, lat0 = None, sizing = 1, \
            left_margin = 0, bottom_margin = 0):
        """
        Sets up the spacial parameters for plotting.
        doc http://gmt.soest.hawaii.edu/doc/5.1.0/gmt.html#j-full
        proj: GMT projection eg 'X' = cartesian, 'M|m' = mercator
        region: tuple containing x_min, x_max, y_min, y_max
        lon0: standard meridian (not always necessary)
        lat0: standard parallel (not always necessary)
        sizing: either scale: distance / degree longitude at meridian
                    or width: total distance of region
        left_margin: move plotting origin in the X direction
        right_margin: move plotting origin in the Y direction
        """
        # work out projection format
        if lon0 == None:
            gmt_proj = '-J%s%s' % (proj, sizing)
        elif lat0 == None:
            gmt_proj = '-J%s%s/%s' % (proj, lon0, sizing)
        else:
            gmt_proj = '-J%s%s/%s/%s' % (proj, lon0, lat0, sizing)

        cmd = [GMT, 'psxy', '-T', gmt_proj, '-X%s' % (left_margin), \
                '-Y%s' % (bottom_margin), '-K', \
                '-R%s/%s/%s/%s' % region]
        # one of the functions that can be run on a blank file
        # as such, '-O' flag needs to be taken care of
        if self.new:
            self.new = False
        else:
            cmd.append('-O')

        Popen(cmd, stdout = self.psf, cwd = self.wd).wait()

    def text(self, x, y, text, dx = 0, dy = 0, align = 'CB', \
            size = '10p', font = 'Helvetica', colour = 'black'):
        """
        Add text to plot.
        x: x position
        y: y position
        text: text to add
        dx: x position offset
        dy: y position offset
        align: Left Centre Right, Top, Middle, Bottom
        size: font size
        font: font familly
        colour: font colour
        """
        tproc = Popen([GMT, 'pstext', '-J', '-R', '-K', '-O', \
                '-D%s/%s' % (dx, dy), '-N', \
                '-F+f%s,%s,%s+j%s' % (size, font, colour, align)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        tproc.communicate('%s %s %s\n' % (x, y, text))
        tproc.wait()

    def sites(self, site_names, shape = 'c', size = 0.1, \
            width = 0.8, colour = 'black', \
            fill = 'gainsboro', transparency = 50, spacing = 0.08, \
            font = 'Helvetica', font_size = '10p', font_colour = 'black'):
        """
        Add sites to map.
        site_names: list of sites to add from defined dictionary
            append ',LB' to change alignment to 'LB' or other
        """
        # step 1: add points on map
        sites_xy = '\n'.join([' '.join(map(str, sites[x.split(',')[0]][:2])) \
                for x in site_names])
        sproc = Popen([GMT, 'psxy', '-J', '-R', '-S%s%s' % (shape, size), \
                '-G%s@%s' % (fill, transparency), '-K', '-O', \
                '-W%s,%s' % (width, colour)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        sproc.communicate(sites_xy)
        sproc.wait()

        # step 2: label points
        # array of x, y, alignment, name
        xyan = []
        for i, xy in enumerate(sites_xy.split('\n')):
            try:
                # user has decided to override position
                name, align = site_names[i].split(',')
            except ValueError:
                # using default position
                name = site_names[i]
                align = sites[name][2]
            xyan.append('%s %s %s' % (xy, align, name))

        tproc = Popen([GMT, 'pstext', '-J', '-R', '-K', '-O', \
                '-Dj%s/%s' % (spacing, spacing), \
                '-F+j+f%s,%s,%s+a0' % (font_size, font, font_colour)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        tproc.communicate('\n'.join(xyan))
        tproc.wait()

    def water(self, colour = 'lightblue', res = 'f'):
        """
        Adds water areas.
        colour: colour of water
        res: resolution
        """
        # GMT land areas are made up of smaller segments
        # as such you can see lines on them and affect visuals
        # therefore the entire area is filled, but then clipped to water
        # pscoast etc can also slightly overlay tickmark (map) outline

        # start cropping to only show wet areas
        Popen([GMT, 'pscoast', '-J', '-R', '-D%s' % (res), \
                '-Sc', '-K', '-O'], \
                stdout = self.psf, cwd = self.wd).wait()
        # fill land and water to prevent segment artifacts
        Popen([GMT, 'pscoast', '-J', '-R', '-G%s' % (colour), \
                '-D%s' % (res), '-K', '-O', '-S%s' % (colour)], \
                stdout = self.psf, cwd = self.wd).wait()
        # crop (-Q) land area off to show only water
        Popen([GMT, 'pscoast', '-J', '-R', '-Q', '-K', '-O'], \
                stdout = self.psf, cwd = self.wd).wait()

    def land(self, fill = 'lightgray', res = 'f'):
        """
        Fills land area.
        fill: colour of land
        res: resolution 'f' full, 'h' high, 'i' intermediate, 'l' low, 'c' crude
        """
        # just like with water, land will show segment artifacts
        # therefore the whole area needs to be filled
        # then cropped to only include land
        Popen([GMT, 'pscoast', '-J', '-R', '-D%s' % (res), '-G%s' % (fill), \
                '-K', '-O'], stdout = self.psf, cwd = self.wd).wait()

    def topo(self, topo_file, topo_file_illu = None, cpt = 'gray'):
        """
        Creates a topography surface using topo files and a colour palette.
        topo_file: file containing topography data
        topo_file: file containing illumination data corresponding to topo_file
            usually the same filename ending with '_i5'
            if not given then the above rule is assumed
        cpt: colour palette to use to display height
        """
        # assume illumination file if not explicitly given
        # assuming the last part of the file is a file extention
        if topo_file_illu == None:
            parts = topo_file.split('.')
            parts[-2] += '_i5'
            topo_file_illu = '.'.join(parts)

        # Q here makes NaN transparent
        Popen([GMT, 'grdimage', topo_file, '-I%s' % (topo_file_illu), '-C%s' % (cpt), \
                '-J', '-R', '-K', '-O', '-Q'], stdout = self.psf, cwd = self.wd).wait()

    def coastlines(self, width = 0.3, colour = 'black', res = 'f'):
        """
        Draws outline of land.
        width: thickness of line
        colour: colour of line
        res: resolution of coastlines
        """
        Popen([GMT, 'pscoast', '-J', '-R', '-D%s' % (res), '-K', '-O', \
                '-W%s,%s' % (width, colour)], \
                stdout = self.psf, cwd = self.wd).wait()

    def ticks(self, major = '60m', minor = '30m', sides = 'ws'):
        """
        Draws map ticks around the edge.
        Note if map doesn't have a left or bottom margin, these will be cut.
        Also part of the map outline may be drawn over by land and/or water.
        It is advisable therefore that ticks are added after area is finished.
        major: these increments have a longer tick
        minor: these increments have a short tick only
        sides: major increments on these sides are labeled with text
        """
        # add sides which aren't wanted as all have to be present
        sides = sides.upper()
        for direction in ['N', 'E', 'S', 'W']:
            if direction not in sides:
                sides = '%s%s' % (sides, direction.lower())

        # here -N3 will also draw marine boundaries
        # they don't show around NZ but something else has to be drawn
        Popen([GMT, 'pscoast', '-J', '-R', '-K', '-O', '-N3', \
                '-Ba%sf%s%s' % (major, minor, sides)], \
                stdout = self.psf, cwd = self.wd).wait()

    def points(self, xy_file, shape = 't', size = 0.08, \
            fill = None, line = 'white', line_thickness = '0.8p'):
        """
        Adds points to map.
        xy_file: file containing x, y positions to plot
        shape: shape to plot at positions
        size: size of shape
        fill: fill colour of shape (default transparent)
        line: line colour of shape
        line_thickness: how thick the outline is
        """
        # build command based on optional fill and thickness
        cmd = [GMT, 'psxy', '-J', '-R', xy_file, \
                '-S%s%s' % (shape, size), '-K', '-O']
        if fill != None:
            cmd.append('-G%s' % (fill))
        if line != None:
            cmd.append('-W%s,%s' % (line_thickness, line))

        Popen(cmd, stdout = self.psf, cwd = self.wd).wait()

    def path(self, in_data, is_file = True, close = False, \
            width = '0.4p', colour = 'black', split = None):
        """
        Draws a path between points.
        in_data: either a filepath to file containing x, y points
                    or a string containing the x, y points
        is_file: whether in_data is a filepath (True) or a string (False)
        close: whether to close the path by joining the first and last points
        width: thickness of line
        colour: colour of line
        split: None continuous, '-' dashes, '.' dots
        """
        # build command based on parameters
        pen = '-W%s,%s' % (width, colour)
        if split != None:
            pen = '%s,%s' % (pen, split)
        cmd = [GMT, 'psxy', '-J', '-R', '-K', '-O', pen]
        if close:
            cmd.append('-L')

        if is_file:
            cmd.append(in_data)
            Popen(cmd, stdout = self.psf, cwd = self.wd).wait()
        else:
            p = Popen(cmd, stdin = PIPE, stdout = self.psf, cwd = self.wd)
            p.communicate(in_data)
            p.wait()

    def seismo(self, src, time, fmt = 'time', \
            width = '1p', colour = 'red'):
        """
        Plots seismograms on map.
        src: file contaning the seismogram data
        time: draw the seismogram up to this reading
        fmt: format of the src file
            'inc' values are read sequentially
            'time' values are read by time
        width: width of the seismo line
        colour: colour of the seismo line
        """
        # grep much faster than python
        # wd same as for GMT for consistency
        if fmt == 'time':
            gp = Popen(['grep', src, '-e', '^>TS%d ' % (time), \
                    '-A%d' % (time + 1)], stdout = PIPE, cwd = self.wd)
        elif gmt == 'inc':
            gp = Popen(['grep', src, '-e', '^>', \
                    '-A%d' % (time + 1)], stdout = PIPE, cwd = self.wd)
        gmt_in = gp.communicate()[0]
        gp.wait()

        sp = Popen([GMT, 'psxy', '-J', '-R', '-N', '-K', '-O',
                '-W%s,%s' % (width, colour)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        sp.communicate(gmt_in)
        sp.wait()

    def cpt_scale(self, x, y, cpt, major, minor, label = None, \
            length = 4.0, thickness = 0.15, horiz = True, \
            arrow_f = True, arrow_b = False):
        """
        Draws a colour palette legend.
        x: x position to place scale
        y: y position to place scale
        cpt: cpt to make scale for
        major: major tick increment (labeled)
        minor: minor tick increment (not labeled)
        label: text label next to scale
        length: how long to draw the scale
        thickness: how thick the scale should be drawn
        horiz: whether to make it horizontal (True) or vertical (False)
        arrow_f: show the forwards continuation arrow (above range)
        arrow_b: show the backwards continuation arrow (belaw range)
        """
        # build command based on parameters
        pos = '-D%s/%s/%s/%s' % (x, y, length, thickness)
        if horiz:
            pos = '%s%s' % (pos, 'h')
        annotation = '-Ba%sf%s' % (major, minor)
        if label != None:
            annotation = '%s:%s:' % (annotation, label.replace(':', ''))
        cmd = [GMT, 'psscale', '-C%s' % (cpt), pos, annotation, '-K', '-O']
        if arrow_f or arrow_b:
            cmd.append('-E%s%s' % ('f' * arrow_f, 'b' * arrow_b))

        Popen(cmd, stdout = self.psf, cwd = self.wd).wait()

    def overlay(self, xyv_file, cpt, dx = '1k', dy = '1k', \
            min_v = None, max_v = None, crop_grd = None, \
            custom_region = None, transparency = 40, climit = 1.0, \
            limit_low = None, limit_high = None, contours = None, \
            contour_thickness = 0.2, contour_colour = 'black', \
            land_crop = False):
        """
        Plot a GMT overlay aka surface.
        xyv_file: file containing x, y and amplitude values
        cpt: cpt to use to visualise data
        dx: x resolution of the surface grid (lower = better quality)
        dy: y resolution of the surface grid
            default unit is longitude/latitude, k: kilometre, e: metre
        min_v: either crop anything below this value (set to NaN)
                   or crop anything above less than max_v if max_v set
        max_v: crop anything below this value that is above min_v
        crop_grd: GMT grd file containing wanted area = 1
        custom_region: grd area region, tuple(x_min, x_max, y_min, y_max)
                speedup is achieved by using a smaller region
        transparency: 0 opaque through 100 invisible
        climit: convergence limit: increasing can drastically improve speed
                if iteration diff is lower than this then result is kept
        limit_low: values below this will be equal to this
        limit_high: values abave this will be equal to this
                limits are one way to make sure values fit in CPT range
                it may be faster to use Numpy pre-processing
        contours: display contour lines every set value or None
        contour_thickness: thickness of contour lines
        contour_colour: colour of contour lines
        land_crop: crop overlay to land area
        """
        # name of intermediate file being worked on
        temp_grd = '%s/%s.grd' % (self.wd, os.path.basename(xyv_file))

        # because we allow setting '-R', backup history file to reset after
        if custom_region != None:
            copyfile(os.path.join(self.wd, 'gmt.history'), \
                    os.path.join(self.wd, 'gmt.history.preserve'))
            region = '-R%s/%s/%s/%s' % custom_region
        else:
            region = '-R'

        # create surface grid
        cmd = [GMT, 'surface', xyv_file, '-G%s' % (temp_grd), '-T0.0', \
                '-I%s/%s' % (dx, dy), '-bi3f', \
                '-C%s' % (climit), region]
        if limit_low != None:
            cmd.append('-Ll%s' % (limit_low))
        if limit_high != None:
            cmd.append('-Lu%s' % (limit_high))
        # ignore stderr: usually because no data in area
        # algorithm in 'surface' is known to fail (no output)
        for attempt in xrange(5):
            Popen(cmd, stderr = self.sink, cwd = self.wd).wait()
            if os.path.exists(temp_grd):
                break
            else:
                print('creating overlay grd attempt %d failed. trying again.' \
                        % (attempt + 1))
        if not os.path.exists(temp_grd):
            print('failed to create grd from %s. no overlay produced.' \
                    % (os.path.basename(xyv_file)))
            return

        # crop to path area by grd file
        if crop_grd != None:
            Popen([GMT, 'grdmath', temp_grd, crop_grd, 'MUL', '=', temp_grd], \
                    cwd = self.wd).wait()

        # crop minimum value
        if min_v != None and max_v == None:
            # ignore stderr: usually because no data in area
            Popen([GMT, 'grdclip', temp_grd, '-G%s' % (temp_grd), \
                    '-Sb%s/NaN' % (min_v)], stderr = self.sink, \
                    cwd = self.wd).wait()

        # restore '-R' if changed
        if custom_region != None:
            move(os.path.join(self.wd, 'gmt.history.preserve'), \
                    os.path.join(self.wd, 'gmt.history'))

        # clip path for land to crop overlay
        if land_crop:
            Popen([GMT, 'pscoast', '-J', '-R', '-Df', '-Gc', \
                    '-K', '-O'], stdout = self.psf, cwd = self.wd).wait()

        # add resulting grid onto map
        # here '-Q' will make NaN transparent
        cmd = [GMT, 'grdimage', temp_grd, '-J', '-R', '-C%s' % (cpt), \
                '-t%s' % (transparency), '-Q', '-K', '-O']
        # ignore stderr: usually because no data in area
        Popen(cmd, stdout = self.psf, stderr = self.sink, \
                cwd = self.wd).wait()

        # add contours
        if contours != None:
            Popen([GMT, 'grdcontour', '-J', '-R', temp_grd, \
            '-C%s' % (contours), '-K', '-O', \
            '-W%s,%s' % (contour_thickness, contour_colour)], \
                    stdout = self.psf, stderr = self.sink, \
                    cwd = self.wd).wait()

        # apply land clip path
        if land_crop:
            Popen([GMT, 'pscoast', '-J', '-R', '-Q', '-K', '-O'], \
                    stdout = self.psf, cwd = self.wd).wait()

        # grd file not needed anymore, prevent clutter
        os.remove(temp_grd)

    def fault(self, in_path, is_srf = False, \
            hyp_shape = 'a', hyp_size = 0.35, \
            plane_width = '1p', plane_colour = 'black', \
            top_width = '2p', top_colour = 'black', \
            hyp_width = '1p', hyp_colour = 'black'):
        """
        Plot SRF fault plane onto map.
        Requires shared_srf.py, replaces addStandardFaultPlane.sh
        in_path: location of input file
        is_srf: if True, input is SRF file. if False, is Corners file.
        hyp_shape: shape to plot at hypocentre 'a' for a star
        hyp_size: size of hypocentre shape
        plane_width: width of line making up fault planes
        plane_colour: colour of line making up fault planes
        top_width: as above for the top edge
        top_colour: as above for the top edge
        hyp_width: as above for hyp_shape outline
        hyp_colour: as above for hyp_shape outline
        """
        if is_srf:
            # use SRF library to retrieve info
            bounds = get_bounds(in_path)
            hypocentre = get_hypo(in_path)

            # process for input into GMT
            gmt_bounds = [['%s %s' % tuple(corner) for corner in plane] \
                    for plane in bounds]
            top_edges = '\n>\n'.join(['\n'.join(corners[:2]) \
                    for corners in gmt_bounds])
            all_edges = '\n>\n'.join(['\n'.join(corners) \
                    for corners in gmt_bounds])
            hypocentre = '%s %s' % tuple(hypocentre)
        else:
            # standard corners file
            bounds = []
            corners = []
            with open(in_path) as cf:
                for line in cf:
                    if line[0] != '>':
                        # not a comment
                        corners.append(line)
                    elif len(corners):
                        # break in long lat stream
                        bounds.append(corners)
                        corners = []
                bounds.append(corners)

            # process for input into GMT
            hypocentre = bounds[0][0]
            top_edges = '>\n'.join([''.join(c[:2]) for c in bounds[1:]])
            all_edges = '>\n'.join([''.join(c) for c in bounds[1:]])

        # plot planes
        planep = Popen([GMT, 'psxy', '-J', '-R', '-L', '-K', '-O', \
                '-W%s,%s,-' % (plane_width, plane_colour)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        planep.communicate(all_edges)
        planep.wait()
        # plot top edges
        topp = Popen([GMT, 'psxy', '-J', '-R', '-K', '-O', \
                '-W%s,%s' % (top_width, top_colour)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        topp.communicate(top_edges)
        topp.wait()
        # hypocentre
        hypp = Popen([GMT, 'psxy', '-J', '-R', '-K', '-O', \
                '-W%s,%s' % (hyp_width, hyp_colour), \
                '-S%s%s' % (hyp_shape, hyp_size)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        hypp.communicate(hypocentre)
        hypp.wait()

    def finalise(self):
        """
        Finalises the postscript.
        """
        # finalisation by running a GMT command without '-K'
        Popen([GMT, 'psxy', '-J', '-R', '-O', '-T'], \
                stdout = self.psf, cwd = self.wd).wait()
        # no more modifications allowed
        self.psf.close()

    def leave(self):
        """
        Alternative to finalise where the file is only closed.
        Useful if this file is opened later.
        """
        self.psf.close()

    def enter(self):
        """
        Only used after leave. Opens file again to continue editing.
        Useful if file is to be externally modified in-between.
        """
        self.psf = open(self.pspath, 'a')

    def png(self, out_dir = None, dpi = 96, clip = True):
        """
        Renders a PNG from the PS.
        Unfortunately relatively slow.
        Could be modified for more formats if needed.
        out_dir: folder to put output in (name as input, different extention)
        dpi: pixels per inch
        clip: whether to crop all whitespace
        """
        # default to output in same directory
        if out_dir == None:
            out_dir = os.path.dirname(self.pspath)

        # A pspath only containing a filename would result in ''
        if out_dir == '':
            out_dir = '.'

        cmd = [GMT, psconvert, self.pspath, '-TG', \
                '-E%s' % (dpi), '-D%s' % (out_dir)]
        if clip:
            cmd.append('-A')
        call(cmd)


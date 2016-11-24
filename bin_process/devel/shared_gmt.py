"""
Library which simplifies the use of GMT.
Functions should hide obvious parameters.
"""

import os
from subprocess import call, PIPE, Popen

# if not in $PATH, can specify location manually
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
    ps2png = 'ps2raster'
else:
    ps2png = 'psconvert'

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

def gmt_defaults(font_annot_primary = 16, map_tick_length_primary = '0.05i', \
        font_label = 16, ps_page_orientation = 'portrait', map_frame_pen = '1p', \
        format_geo_map = 'D', map_frame_type = 'plain', format_float_out = '%lg', \
        proj_length_unit = 'i', ps_media = 'A2', extra = []):
    """
    Sets default values for GMT.
    GMT stores these values in the file 'gmt.defaults'
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
    call(cmd)

###
### MAIN PLOTTING CLASS
###
class GMTPlot:

    def __init__(self, pspath):
        self.pspath = pspath
        self.psf = open(pspath, 'w')
        # figure out where to run GMT from
        self.wd = os.path.dirname(pspath)
        if self.wd == '':
            self.wd = '.'
        # if previous/custom gmt.defaults and gmt.history was wanted
        # it should have already been copied into this directory
        if not os.path.exists(os.path.join(self.wd, 'gmt.defaults')):
            gmt_defaults()

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
        proc = Popen([GMT, 'psxy', '-K', '-G%s' % (colour), \
                '-JX%s/%s' % (length, height), '-R0/%s/0/%s' % (length, height), \
                '-Xa-%s' % (left_margin), '-Ya-%s' % (bottom_margin)], \
                stdin = PIPE, stdout = self.psf, cwd = self.wd)
        proc.communicate('%s 0\n%s %s\n0 %s\n0 0' \
                % (length, length, height, height))
        proc.wait()

    def spacial(self, proj, x_min, x_max, y_min, y_max, \
            lon0 = None, lat0 = None, sizing = 1, \
            left_margin = 0, bottom_margin = 0):
        """
        Sets up the spacial parameters for plotting.
        doc http://gmt.soest.hawaii.edu/doc/5.1.0/gmt.html#j-full
        proj: GMT projection eg 'X' = cartesian, 'M|m' = mercator
        x_min: minimal X position
        x_max: maximal X position
        y_min: minimal Y position
        y_max: maximal Y position
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

        Popen([GMT, 'psxy', '-T', gmt_proj, '-X%s' % (left_margin), '-Y%s' % (bottom_margin), \
                '-O', '-K', '-R%s/%s/%s/%s' % (x_min, x_max, y_min, y_max)], \
                stdout = self.psf, cwd = self.wd).wait()

    def land(self, fill = 'lightgray', res = 'f'):
        """
        Fills land area.
        fill: colour of land
        res: resolution 'f' full, 'h' high, 'i' intermediate, 'l' low, 'c' crude
        """
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
                '-W%s,%s' % (width, colour)], stdout = self.psf, cwd = self.wd).wait()

    def water(self, colour = 'lightblue', res = 'f'):
        """
        Adds water areas.
        colour: colour of water
        res: resolution
        """
        Popen([GMT, 'pscoast', '-J', '-R', '-D%s' % (res), '-K', '-O', \
                '-S%s' % (colour)], stdout = self.psf, cwd = self.wd).wait()

    def ticks(self, major = '60m', minor = '30m', sides = 'ws'):
        """
        Draws map ticks around the edge.
        Note if map doesn't have a left or bottom margin, these will be cut.
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

    def path(self, indata, isfile = True, close = False, \
            width = '0.4p', colour = 'black', split = None):
        """
        Draws a path between points.
        indata: either a filepath to file containing x, y points
                    or a string containing the x, y points
        isfile: whether indata is a filepath (True) or a string (False)
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

        if isfile:
            cmd.append(indata)
            Popen(cmd, stdout = self.psf, cwd = self.wd).wait()
        else:
            p = Popen(cmd, stdin = PIPE, stdout = self.psf, cwd = self.wd)
            p.communicate(indata)
            p.wait()

    def cpt_scale(self, x, y, cpt, major, minor, label = None, \
            length = 4.0, thickness = 0.15, horiz = True):
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
        """
        # build command based on parameters
        pos = '-D%s/%s/%s/%s' % (x, y, length, thickness)
        if horiz:
            pos = '%s%s' % (pos, 'h')
        annotation = '-Ba%sf%s' % (major, minor)
        if label != None:
            annotation = '%s:%s:' % (annotation, label.replace(':', ''))
        cmd = [GMT, 'psscale', '-C%s' % (cpt), pos, annotation, '-K', '-O']

        Popen(cmd, stdout = self.psf, cwd = self.wd).wait()

    def finalise(self):
        """
        Finalises the postscript.
        """
        # finalisation by running a GMT command without '-K'
        Popen([GMT, 'psxy', '-J', '-R', '-O', '-T'], \
                stdout = self.psf, cwd = self.wd).wait()
        # no more modifications allowed
        self.psf.close()

    def png(self, out_dir = None, resolution = 120, clip = True):
        """
        Renders a PNG from the PS.
        Could be modified for more formats if needed.
        out_dir: folder to put output in (name as input, different extention)
        resolution: DPI
        clip: whether to crop all white
        """
        # default to output in same directory
        if out_dir == None:
            out_dir = os.path.dirname(self.pspath)

        # A pspath only containing a filename would result in ''
        if out_dir == '':
            out_dir = '.'

        cmd = [GMT, ps2png, self.pspath, '-TG', \
                '-E%s' % (resolution), '-D%s' % (out_dir)]
        if clip:
            cmd.append('-A')
        call(cmd)


#!/usr/bin/env python

import os
import sys
from glob import glob
from shutil import copyfile
from subprocess import call, Popen, PIPE
import os.path
sys.path.append(os.path.abspath(os.path.curdir))

# make temp file with the model coord boundaries
mp_tmp = 'tmp.modelpath'
if os.path.isfile(mp_tmp):
    os.remove(mp_tmp)
if os.path.exists(mp_tmp):
    raise IOError('File `tmp.modelpath` cannot be removed!');
modparms_handle = open(MODELPARAMS, 'r')
mp_handle = open(mp_tmp, 'w')
# transfer corner coordinates
corner_defs = []
for line in modparms_handle:
    if any(corner in line for corner in ['c1=', 'c2=', 'c3=', 'c4=']):
        corner_defs.append(filter(None, line.split(' '))) # eg. ['c1=', 'lon_value', 'lat_value']
for corner in sorted(corner_defs):
    mp_handle.write(' '.join(corner[1:]))
mp_handle.close()
modparms_handle.close()

# create masks for later plotting of different layers
call(['gmt', 'grdlandmask', '-R' + TS_REGION, '-Dh', '-I' + DX + '/' + DY, '-Glandmask.grd'])
call(['gmt', 'grdmask', mp_tmp, '-R' + TS_REGION, '-I' + DX + '/' + DY, '-Gmodelmask.grd'])
call(['gmt', 'grdmath', 'landmask.grd', 'modelmask.grd', 'MUL = allmask.grd'])

# create color palette for plotting the topography
BASECPT = 'y2r_brown.cpt'
TOPOMAX = '4000'
AINC = '10'
AMAX = '50'
AMIN = '2' # min value to show
ABELOW = 'NaN' # value to show when less than AMIN

# color palette for velocity
cpt_handle = open(BASECPT, 'w')
call(['gmt', 'makecpt', '-Chot', '-I', '-T0/' + AMAX + '/' + AINC, '-A50'], stdout = cpt_handle)
cpt_handle.close()

# end of creating color palette

# get the avg Lon and Lat values (i.e. the midpoint of the TSlice image for the geo projection)
AVGLL = ['%0.6f'%((float(PLOT_XMIN) + float(PLOT_XMAX)) / 2), \
        '%0.6f'%((float(PLOT_YMIN) + float(PLOT_YMAX)) / 2)]
# store plot projection and region in GMT format for ease of use later
ATT = '-JT' + '/'.join[AVGLL[0], AVGLL[1], XINCH] + '-R' + PLOT_REGION

# specify the specific strong motion stations of interest for plotting
STATFILE = os.path.join(MAINDIR, StationInfo, cantstations.ll)

# set all plotting defaults to use
call(['gmt', 'gmtset', 'FONT_ANNOT_PRIMARY 16', 'MAP_TICK_LENGTH_PRIMARY 0.05i', \
        'FONT_LABEL 16', 'PS_PAGE_ORIENTATION PORTRAIT', 'MAP_FRAME_PEN 1p', \
        'FORMAT_GEO_MAP D', 'MAP_FRAME_TYPE plain', 'FORMAT_FLOAT_OUT %lg', \
        'PROJ_LENGTH_UNIT i'])

# start loop for all time slices
tsfcnt = TS_START
while tsfcnt < TS_TOTAL:
    # current TS index
    tsspot = tsfcnt * TS_INC
    # get current time in seconds (first get index in floating point fmt)
    tsfnum = `echo $tsfcnt | gawk '{printf "%.4d\n",$1;}'`
    tt = `echo $tsspot $ORIG_DT | gawk '{printf "%.2f\n",$1*$2;}'`
    # specify the names of the .ps and .png files that will be created
    PLOTFILE = 'PlotFiles/ts-str' + str(tsfnum) + '.ps'
    PNGFILE = PNGDIR + '/ts-str' + str(tsfnum) + '.png'
    # remove any existing files
    if os.path.exists(PLOTFILE):
        os.remove(PLOTFILE)
    if os.path.exists(PNGFILE):
        os.remove(PNGFILE)

    # specify plot and panel size (defaults 8.5 x 11)
    plot_handle = open(PLOTFILE, 'w')
    gmt_pipe = Popen(['gmt', 'psxy', '-JX8.5/11', '-R0/8.5/0/11',  '-L', '-G180/180/180', \
            '-X0',  '-Y0',  '-K'], stdout = plot_handle, stdin = PIPE)
    gmt_pipe.communicate(input='#-W0/180/180/180\n' \
            + '0.3 1.0\n' \
            + '0.3 7.8\n' \
            + '6.5 7.8\n' \
            + '6.5 1.0\n'

    # set the color scale
    call(['gmt psscale', '-C' + BASECPT, '-Ef', '-D3.0/2.0/2.5/0.15h', '-K', '-O', \
            '-Ba' + AINC + 'f' + AINC + ':"ground velocity (cm/s)":'], stdout = plot_handle)

    # specify the X and Y offsets for plotting (I dont really understand this yet)
    gmt_pipe = Popen(['gmt', 'psxy', '-V ', ATT, '-L', '-K', '-O', '-X' + XORG, \
            '-Y' + YORG], stdout = plot_handle, stdin = PIPE)
    gmt_pipe.communicate(input='#-W5/255/255/0\n')

    if option == 1:
        # use 'ts2xyz.exe to get TSlice outout in xyz format
        call([os.path.join(WCC_PROGDIR, 'ts2xyz'), 'infile=' + TSFILE, 'outfile=outf', \ 
                'swap_bytes=' + SWAP_BYTES, 'gridfile=' + GRIDFILE, 'xyts=1', 'scale=' + SCALE, \
                'ts=' + tsspot, 'trv=0', 'dxts=' + DXTS, 'dyts=' + DYTS, 'dzts=1', 'absmax=1' \
                'read_header=1','outbin=1', 'lonlat=' + LONLAT_OUT, 'geoproj=1'])

    elif option == 2:
        TSFILEPREFIX = os.path.join(TS_OUTFILE_DIR, NAME)
        outf = TSFILEPREFIX + '_ts' + str(tsfcnt).zfill(4) + '\n'

    # currently components are not looped over, as using ABSMAX=1
    # loop over the different components
    for comp in COMPS
        if option == 2:
            if SWAP_BYTES == '1':
                # get file in correct format - BB added
                gmt xyz2grd ${outf}.${comp} -Soutf.${comp} -V -Zf
            elif SWAP_BYTES == '0':
                copyfile(outf + '.' + comp, 'outf.' + comp)

        # try a different version of plotting
        # clippath for land
        call(['gmt', 'pscoast', ATT, '-Df', '-Gc', '-K', '-O'], stdout = plot_handle)
        # land
        call(['gmt', 'grdimage', TOPO_FILE, ILLU, PALETTE, ATT, '-K' '-O'], stdout = plot_handle)
        # clear clippath
        call(['gmt', 'pscoast', '-R', '-J', '-O', '-K', '-Q'], stdout = plot_handle)

        # add urban areas
        URBANDIR = MAINDIR + '/PlottingData/sourcesAndStrongMotionStations'
        call(['gmt', 'psxy', URBANDIR + '/ChchUrbanBoundary.xy', 'ATT', '-G160/160/160', '-W0.5p', \
                '-O', '-K'], stdout = plot_handle)

        # create ground motion intensity surface from the TSlice output
        call(['gmt', 'surface', 'outf.' + comp, '-Gtmp0.grd', '-I' + DX + '/' + DY, \
                '-R' + TS_REGION, '-T0.0', '-bi3f']) # make surface
        call(['gmt', 'grdclip', 'tmp0.grd', '-Gtmp1.grd', '-Sb' + AMIN + '/' + ABELOW]) # clip minimum
        call(['gmt', 'grdmath', 'modelmask.grd', 'tmp1.grd', \
                'MUL = outf_' + comp + '.grd']) # clip to TS region
        call(['gmt', 'grdclip', 'tmp1.grd', '-Goutf_' + comp + '.grd', \
                '-Sb' + AMIN + '/' + ABELOW]) # clip minimum

        copyfile('outf_' + comp + '.grd', 'tmp1.grd')
        call(['gmt', 'grdmath', 'modelmask.grd', '1', 'SUB', 'tmp1.grd', \
                'ADD = outf_' + comp + '.grd'])
        # add grid image to ps plot
        call(['gmt', 'grdimage', 'outf_' + comp + '.grd',  ATT, '-C' + BASECPT, '-Q', '-t50', \
                '-K', '-O'], stdout = plot_handle)
        # add coastline
        call(['gmt', 'pscoast', '-A0/0/1', '-N1', '-N2', ATT, '-Df', '-S135/205/250', '-W1,black', \
                '-K', '-O'], stdout = plot_handle)
        call(['gmt', 'pscoast', '-A0/2/2', ATT, '-Df', '-W1,black', '-K', '-O'], stdout = plot_handle)

        # add fault planes
        ADDFLTPLANEDIR = MAINDIR + '/PlottingData/sourcesAndStrongMotionStations'
        FAULTFILE = ADDFLTPLANEDIR + '/22Feb2011FaultPlane.xy;'
        LINE = '-W0.5p,black,-'
        TOPEDGE = '-W2p,black'
        HYPOPEN = '-W1p,black;'
        # call fault plane routine
        call(['bash', ADDFLTPLANEDIR + '/addStandardFaultPlane.sh', PLOTFILE, '-R' + TS_REGION, \
                '-JT' + AVGLL[0] + '/' + AVGLL[1] + '/' + XINCH, FAULTFILE, LINE, TOPEDGE, HYPOPEN)

        # add Main and subtitles
        gmt_pipe = Popen(['gmt', 'pstext', ATT, '-N', '-O', '-K', '-D0.0/0.35', \
                '-F+f20p,Helvetica-Bold,black+jLB+a0'], stdout = plot_handle, stdin = PIPE)
        gmt_pipe.communicate(input = PLOT_XMIN + ' ' + PLOT_YMAX + ' ' + MAIN_TITLE + '\n')

        # subtitle
        gmt_pipe = Popen(['gmt', 'pstext', ATT, '-N', '-O', '-K', '-D0.0/0.1', \
                '-F+f+j+a0,'], stdout = plot_handle, stdin = PIPE)
        gmt_pipe.communicate(input = PLOT_XMIN + ' ' + PLOT_YMAX \
                + '14,Helvetica,black LB ' + SUB_TITLE + '\n' \
                + PLOT_XMAX + ' ' + PLOT_YMAX + ' 16,Helvetica,black RB t=' + tt + ' sec\n'

        #add scale to show distance
        gmt psbasemap $ATT -L172.50/-43.90/$AVGLL[1]/25.0 -Ba30mf30mWSen -K -O >> $PLOTFILE

        #loop over all the sites/locations to be displayed
        set x = 0
        foreach site ( $SLON )
            @ x ++

            #plot the location as a point
            gmt psxy $ATT -S$SSYM -G$SFIL -W$SLIN -O -K << END >> $PLOTFILE
$SLON[$x] $SLAT[$x]
END

            #add the location name
            gmt pstext $ATT -N -O -K -Dj0.05/0.05 -F+j+f12,Helvetica,black+a0 << END >>  $PLOTFILE
$SLON[$x] $SLAT[$x] $SPOS[$x] $SITES[$x]
END
        end

        #plot all the strong motion station locations
        gmt psxy $STATFILE $ATT -St0.08 -G000/000/000 -W$SLIN -O -K >> $PLOTFILE
        #shift plotting origin (for 3 component plotting)
        gmt psxy -V $ATT -L -W5,255/255/0 -O -K -X$XSHFT << END >>  $PLOTFILE
END
        #end of component
    end

    #finalize the plot (i.e. no -K)
    gmt psxy -V $ATT -L -W5,255/255/0 -O << END >>  $PLOTFILE
END
    #convert .ps file to raster (.png) file and remove plotfile
    gmt ps2raster $PLOTFILE -A -TG -E$RES -D$PNGDIR
#    \rm $PLOTFILE

    #increment the counter for the TSlice index
    tsfcnt += 1
end

# remove all temporary files
for file in ['tmp0.grd', 'tmp1.grd', 'modelmask.grd', 'landmask.grd', 'allmask.grd', \
        'tmp.modelpath', BASECPT, 'gmt.conf', 'gmt.history'] + glob('outf.?') + glob('outf_?.grd'):
    os.remove(file)


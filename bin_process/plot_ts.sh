#!/usr/bin/env sh

# following variables are used from e3d.par:
#   extended_name, plot_main_title, plot_sub_title, plot_option, wcc_prog_dir,
#   stat_file, vel_mod_params_dir, h, dt_ts, dx_ts, dy_ts, dz_ts, ts_start, ts_inc,
#   swap_bytes, lonlat_out, scale
source e3d.par

#High-level directories
MAINDIR=/hpc/scratch/nesi00213

#information required for option 1.
#----------------------------------
#location of the TSlice binary file created by 'merge_tsP3.csh'
TSFILE=./OutBin/${extended_name}_xyts.e3d
#

#information required for option 2.
#----------------------------------
#specify the location of the TS files
TSFILEDIR=./TSFiles

#details of the spatial and termporal discretization and spacing
ORIG_DT=0.1   #time step of the time slice output (this is DT*DT_TS from the run files)

# used in 'ts2xyz' bin, may not be used
GRIDFILE="${vel_mod_params_dir}/gridout_nz01-h${h}"   #used in 'ts2xyz' exe - so not presentely used
MODELPARAMS="${vel_mod_params_dir}/model_params_nz01-h${h}" #used to get location of model corners

#components to consider (0= when only looking at vector magnitude [i.e. ABSMAX=1 further down].
COMPS=( 0 )

#specify the directory for PNG files and the resolution (dpi) to create
PSDIR=PlotFiles
PNGDIR=Png
RES=140   #720 for PDF quality
#make directories for storing the .ps and .png files
mkdir -p $PSDIR $PNGDIR

TOPODIR=${MAINDIR}/PlottingData/TopoData
TOPO_FILE=${TOPODIR}/srtm_71_21.grd
ILLU=-I${TOPODIR}/srtm_71_21_i5.grd
TOPO_FILE2=${TOPODIR}/etopo2.grd
#set PALETTE=-Crelief.cpt
PALETTE=-Cgray

#specific locations to display on figure (not currently used)
SITES=( "Rolleston" "Darfield" "Lyttelton" "Akaroa" "Kaiapoi" "Rakaia" "Oxford" )
#set SPOS = ( 1 1 1 1 1 1 1)
SPOS=( "RB" "CB" "LM" "RB" "LB" "RT" "LB" )  #alignment; 2letter, L,C,R (left, center, right); T,M,B (top, middle, bottom)
SLON=( 172.3791667 172.1116667 172.7194444 172.9683333 172.6569444 172.0230556 172.1938889 )
SLAT=( -43.59083333 -43.48972222 -43.60305556 -43.80361111 -43.38277778 -43.75611111 -43.29555556 )
#specifying plotting preferences for site locations
SSYM=c0.10 #symbol
SFIL=220/220/220 #fill color
SLIN=1,000/000/000 #line?

#location of offset plotting (for when 3 component plotting used - not currently utilized)
XORG=1.15
YORG=2.5
XINCH=5.0
XSHFT=`echo $XINCH 0.2 | gawk '{printf "%f\n",$1+$2;}'`

#specify the maximum plotting window for the basemap
PLOT_XMIN=171.75
PLOT_XMAX=173.00
PLOT_YMIN=-44.00
PLOT_YMAX=-43.20
#region in GMT format
PLOT_REGION="${PLOT_XMIN}/${PLOT_XMAX}/${PLOT_YMIN}/${PLOT_YMAX}"

#specify the plotting region for the time slice - this must 'within' the time slice output range (check by looking at MODELPARAMS file)
TS_XMIN=171.75
TS_XMAX=173.00
TS_YMIN=-44.00
TS_YMAX=-43.20
#region in GMT format
TS_REGION="${TS_XMIN}/${TS_XMAX}/${TS_YMIN}/${TS_YMAX}"
#specify the increments of X/Y (cartesian coords) for masks etc.
DX=0.002
DY=0.002

#make temp file with the model coord boundaries
\rm tmp.modelpath
foreach corner ( c1 c2 c3 c4 )
    grep "$corner= " $MODELPARAMS | gawk ' { print $2, $3 } ' >> tmp.modelpath
end

#create masks for later plotting of different layers
gmt grdlandmask -R$TS_REGION -Dh -I$DX/$DY -Glandmask.grd
gmt grdmask tmp.modelpath -R$TS_REGION -I$DX/$DY -Gmodelmask.grd
gmt grdmath landmask.grd modelmask.grd MUL = allmask.grd

#create color palette for plotting the topography
BASECPT = y2r_brown.cpt
TOPOMAX = 4000
AINC = 10
AMAX = 80
AMIN = 2 #the min value to show
ABELOW = NaN #the value to show when less than AMIN

# color palette for velocity
gmt makecpt -Chot -I -T0/$AMAX/$AINC -A50 > $BASECPT

#end of creating color palette

#get the avg Lon and Lat values (i.e. the midpoint of the TSlice image for the geo projection)
AVGLL = `echo $PLOT_XMIN $PLOT_XMAX $PLOT_YMIN $PLOT_YMAX | gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`
#store plot projection and region in GMT format for ease of use later
ATT = "-JT$AVGLL[1]/$AVGLL[2]/${XINCH} -R${PLOT_REGION}"

#set all plotting defaults to use
gmt gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i

#start loop for all time slices
@ tsfcnt = $ts_start
while ( $tsfcnt < $ts_total )
    #current TS index
    @ tsspot = ( $tsfcnt * $TS_INC )
    #get current time in seconds (first get index in floating point fmt)
    tsfnum = `echo $tsfcnt | gawk '{printf "%.4d\n",$1;}'`
    tt = `echo $tsspot $ORIG_DT | gawk '{printf "%.2f\n",$1*$2;}'`
    #specify the names of the .ps and .png files that will be created
    PLOTFILE = PlotFiles/ts-str${tsfnum}.ps
    PNGFILE = $PNGDIR/ts-str${tsfnum}.png
    #remove any existing files
    \rm $PLOTFILE $PNGFILE

    #specify plot and panel size (defaults 8.5 x 11)
    EDGE_COLOR = 255/255/255 #180/180/180 = grey ; 255/255/255=white
    #gmt psxy -JX8.5/11 -R0/8.5/0/11 -L -G180/180/180 -X0 -Y0 -K << END > $PLOTFILE #-W0/180/180/180
    gmt psxy -JX8.5/11 -R0/8.5/0/11 -L -G${EDGE_COLOR} -X0 -Y0 -K << END > $PLOTFILE #-W0/180/180/180
0.3 1.0
0.3 7.8
6.5 7.8
6.5 1.0
END

    #set the color scale
    gmt psscale -C$BASECPT -Ef -D3.0/2.0/2.5/0.15h -K -O -Ba${AINC}f${AINC}:"ground velocity (cm/s)": >> $PLOTFILE

    #specify the X and Y offsets for plotting (I dont really understand this yet)
    gmt psxy -V $ATT -L  -K -O -X$XORG -Y$YORG << END >> $PLOTFILE #-W5/255/255/0
END

    if [ "$plot_option" -eq 1 ]; then
#    use 'ts2xyz.exe. to get TSlice outout in xyz format
        ${wcc_prog_dir}/ts2xyz infile=$TSFILE outfile=outf swap_bytes=$swap_bytes \
       gridfile=$GRIDFILE xyts=1 scale=$scale ts=$tsspot trv=0 \
       dxts=$dxts dyts=$dyts dzts=1 absmax=1 \
       read_header=1 outbin=1 lonlat=$lonlat_out geoproj=1

    elif [ "$plot_option" -eq 2 ]; then
        TSFILEPREFIX = ${TSFILEDIR}/${extended_name}
        outf = `echo $TSFILEPREFIX $tsfcnt | gawk '{printf "%s_ts%.4d\n",$1,$2;}'`
    fi

#currently components are not looped over, as using ABSMAX=1
    #loop over the different components
    c = 0
    foreach comp ( $COMPS )
        @ c ++
        if [ "$plot_option" -eq 2 ]; then
            if ($swap_bytes == 1) then
                #get file in correct format - BB added
                gmt xyz2grd ${outf}.${comp} -Soutf.${comp} -V -Zf
            else if ($swap_bytes == 0) then
                \cp ${outf}.${comp} outf.${comp}
            endif
        fi

        #try a different version of plotting
        # clippath for land
        gmt pscoast $ATT -Df -Gc -K -O >> $PLOTFILE
        # land
        gmt grdimage $TOPO_FILE $ILLU $PALETTE $ATT -K -O >> $PLOTFILE
        # clear clippath
        gmt pscoast -R -J -O -K -Q >> $PLOTFILE

        #add urban areas
        URBANDIR = ${MAINDIR}/PlottingData/sourcesAndStrongMotionStations
        gmt psxy ${URBANDIR}/ChchUrbanBoundary.xy $ATT -G160/160/160 -W0.5p -O -K >> $PLOTFILE

        #create ground motion intensity surface from the TSlice output
        gmt surface outf.${comp} -Gtmp0.grd -I$DX/$DY -R$TS_REGION -T0.0 -bi3f  #make surface
        gmt grdclip tmp0.grd -Gtmp1.grd -Sb${AMIN}/$ABELOW #clip minimum
        gmt grdmath modelmask.grd tmp1.grd MUL = outf_${comp}.grd #clip to TS region
        gmt grdclip tmp1.grd -Goutf_${comp}.grd -Sb${AMIN}/$ABELOW #clip minimum

        \cp outf_${comp}.grd tmp1.grd
        gmt grdmath modelmask.grd 1 SUB tmp1.grd ADD = outf_${comp}.grd
        #add grid image to ps plot
        gmt grdimage outf_${comp}.grd $ATT -C$BASECPT -Q -t50 -K -O >> $PLOTFILE
        #add coastline
        gmt pscoast -A0/0/1 -N1 -N2 $ATT -Df -S135/205/250 -W1,black -K -O >> $PLOTFILE
        gmt pscoast -A0/2/2 $ATT -Df -W1,black -K -O >> $PLOTFILE

        #Add fault planes
        ADDFLTPLANEDIR = ${MAINDIR}/PlottingData/sourcesAndStrongMotionStations
        FAULTFILE = ${ADDFLTPLANEDIR}/bev01_DarfieldFaultPlane.xy;
        LINE = -W0.5p,black,-
        TOPEDGE = -W2p,black
        HYPOPEN = -W1p,black;
        #call fault plane routine
        bash ${ADDFLTPLANEDIR}/addStandardFaultPlane.sh $PLOTFILE -R$TS_REGION -JT$AVGLL[1]/$AVGLL[2]/${XINCH} $FAULTFILE $LINE $TOPEDGE $HYPOPEN

        #add Main and subtitles
        gmt pstext $ATT -N -O -K -D0.0/0.35 -F+f20p,Helvetica-Bold,black+jLB+a0 << END >>  $PLOTFILE #
$PLOT_XMIN $PLOT_YMAX $plot_main_title
END

        #subtitle
        gmt pstext $ATT -N -O -K -D0.0/0.1 -F+f+j+a0, << END >>  $PLOTFILE
$PLOT_XMIN $PLOT_YMAX 14,Helvetica,black LB $plot_sub_title
$PLOT_XMAX $PLOT_YMAX 16,Helvetica,black RB t=$tt sec
END

        #add scale to show distance
        gmt psbasemap $ATT -L172.50/-43.90/$AVGLL[2]/25.0 -Ba30mf30mWSen -K -O >> $PLOTFILE

        #loop over all the sites/locations to be displayed
        x = 0
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
        gmt psxy "$stat_file" $ATT -St0.08 -G000/000/000 -W$SLIN -O -K >> $PLOTFILE
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
    @ tsfcnt ++
end

#remove all temporary files
\rm tmp0.grd tmp1.grd modelmask.grd landmask.grd allmask.grd
\rm tmp.modelpath outf.? outf_?.grd
\rm $BASECPT
\rm gmt.conf gmt.history
#\rm -r $PLOTFILE

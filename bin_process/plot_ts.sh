#!/usr/bin/env bash
# bash compatible, should work with ksh too, find alternative to arrays for sh compatibility

# following variables are used from e3d.par:
#   extended_name, plot_main_title, plot_sub_title, plot_option, wcc_prog_dir,
#   stat_file, vel_mod_params_dir, h, dt_ts, dx_ts, dy_ts, dz_ts, ts_start, ts_inc,
#   swap_bytes, lonlat_out, scale, plot_sites, plot_s_pos, plot_s_lon, plot_s_lat,
#   plot_s_sym, plot_s_fil, plot_s_lin, plot_x_org, plot_y_org, plot_x_inch, plot_x_shift,
#   plot_x_min, plot_x_max, plot_y_min, plot_y_max, plot_region, plot_ts_x_min, plot_ts_x_max,
#   plot_ts_y_min, plot_ts_y_max, plot_region, plot_ts_region, plot_dx, plot_dy, plot_palette,
#   ts_file, ts_out_prefix, plot_ps_dir, plot_png_dir, plot_res, plot_orig_dt, plot_comps,
#   global_root
source e3d.par

# used by 'ts2xyz' bin, may not be used
grid_file="${vel_mod_params_dir}/gridout_nz01-h${h}"
model_params="${vel_mod_params_dir}/model_params_nz01-h${h}" # to get location corners

# output dirs
mkdir -p $plot_ps_dir $plot_png_dir

TOPODIR=${global_root}/PlottingData/TopoData
TOPO_FILE=${TOPODIR}/srtm_71_21.grd
ILLU=-I${TOPODIR}/srtm_71_21_i5.grd
TOPO_FILE2=${TOPODIR}/etopo2.grd

# make temp file with the model coord boundaries
\rm tmp.modelpath
for corner in c1 c2 c3 c4; do
    grep "$corner= " $model_params | gawk ' { print $2, $3 } ' >> tmp.modelpath
done

# create masks for later plotting of different layers
gmt grdlandmask -R$plot_ts_region -Dh -I$plot_dx/$plot_dy -Glandmask.grd
gmt grdmask tmp.modelpath -R$plot_ts_region -I$plot_dx/$plot_dy -Gmodelmask.grd
gmt grdmath landmask.grd modelmask.grd MUL = allmask.grd

# create color palette for plotting the topography
BASECPT=y2r_brown.cpt
TOPOMAX=4000
AINC=10
AMAX=80
AMIN=2 #the min value to show
ABELOW=NaN #the value to show when less than AMIN

# fault plane
ADDFLTPLANEDIR=${global_root}/PlottingData/sourcesAndStrongMotionStations
FAULTFILE=${ADDFLTPLANEDIR}/bev01_DarfieldFaultPlane.xy;
LINE=-W0.5p,black,-
TOPEDGE=-W2p,black
HYPOPEN=-W1p,black;

# color palette for velocity
gmt makecpt -Chot -I -T0/$AMAX/$AINC -A50 > $BASECPT

#end of creating color palette

# avg Lon/Lat (midpoint of the TSlice image for the geo projection)
avg_ll=(`echo $plot_x_min $plot_x_max $plot_y_min $plot_y_max | gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`)
# plot projection and region in GMT format for ease of use later
att="-JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} -R${plot_region}"


add_site() {
    # plot location as a point
    gmt psxy $att -S$plot_s_sym -G$plot_s_fil -W$plot_s_lin -O -K << END >> $plot_file
${plot_s_lon[$1]} ${plot_s_lat[$1]}
END

    # add location name
    gmt pstext $att -N -O -K -Dj0.05/0.05 -F+j+f12,Helvetica,black+a0 << END >>  $plot_file
${plot_s_lon[$1]} ${plot_s_lat[$1]} ${plot_s_pos[$1]} ${plot_sites[$1]}
END
}


###################### BEGIN TEMPLATE ##########################
echo Creating PS Template...
plot_file_template=plot_template.ps
# specify plot and panel size (defaults 8.5 x 11)
edge_colour=255/255/255 #180/180/180 = grey ; 255/255/255=white
gmt psxy -JX8.5/11 -R0/8.5/0/11 -L -G${edge_colour} -X0 -Y0 -K << END > $plot_file_template #-W0/180/180/180
0.3 1.0
0.3 7.8
6.5 7.8
6.5 1.0
END
# set the color scale
gmt psscale -C$BASECPT -Ef -D3.0/2.0/2.5/0.15h -K -O -Ba${AINC}f${AINC}:"ground velocity (cm/s)": >> $plot_file_template
# specify the X and Y offsets for plotting (I dont really understand this yet)
gmt psxy -V $att -L  -K -O -X$plot_x_org -Y$plot_y_org << END >> $plot_file_template 2>/dev/null #-W5/255/255/0
END
# try a different version of plotting
# clippath for land
gmt pscoast $att -Df -Gc -K -O >> $plot_file_template
# land
gmt grdimage $TOPO_FILE $ILLU $plot_palette $att -K -O >> $plot_file_template
# clear clippath
gmt pscoast -R -J -O -K -Q >> $plot_file_template
# add urban areas
URBANDIR=${global_root}/PlottingData/sourcesAndStrongMotionStations
gmt psxy ${URBANDIR}/ChchUrbanBoundary.xy $att -G160/160/160 -W0.5p -O -K >> $plot_file_template
echo Template Complete.
####################### END TEMPLATE ###########################


# set all plotting defaults to use
gmt gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i

#start loop for all time slices
tsfcnt=$ts_start
while [ $tsfcnt -lt $ts_total ]; do
    # current TS index
    tsspot=$(($tsfcnt * $ts_inc))
    echo Rendering timeslice $tsfcnt/$ts_total \($tsspot\)

    # output number (4 decimal places)
    tsfnum=`echo $tsfcnt | gawk '{printf "%.4d\n",$1;}'`
    # outputs
    plot_file=$plot_ps_dir/ts-str${tsfnum}.ps
    png_file=$plot_png_dir/ts-str${tsfnum}.png
    #remove any existing files
    if [ -e "$png_file" ]; then
        \rm "$png_file"
    fi
    cat "$plot_file_template" > "$plot_file"

    if [ "$plot_option" -eq 1 ]; then
       # use 'ts2xyz.exe. to get TSlice outout in xyz format
        ${wcc_prog_dir}/ts2xyz infile=$ts_file outfile=outf swap_bytes=$swap_bytes \
                gridfile=$grid_file xyts=1 scale=$scale ts=$tsspot trv=0 \
                dxts=$dxts dyts=$dyts dzts=1 absmax=1 \
                read_header=1 outbin=1 lonlat=$lonlat_out geoproj=1
    elif [ "$plot_option" -eq 2 ]; then
        outf=`echo $ts_out_prefix $tsfcnt | gawk '{printf "%s_ts%.4d\n",$1,$2;}'`
    fi

    # currently components are not looped over, as using ABSMAX=1
    # loop over the different components
    for comp in "${plot_comps[@]}"; do
        if [ "$plot_option" -eq 2 ]; then
            if [ "$swap_bytes" -eq 1 ]; then
                # get file in correct format - BB added
                gmt xyz2grd ${outf}.${comp} -Soutf.${comp} -V -Zf
            elif [ "$swap_bytes" -eq 0 ]; then
                \cp ${outf}.${comp} outf.${comp}
            fi
        fi

        #create ground motion intensity surface from the TSlice output
        gmt surface outf.${comp} -Gtmp0.grd -I$plot_dx/$plot_dy -R$plot_ts_region -T0.0 -bi3f 2>/dev/null
        gmt grdclip tmp0.grd -Gtmp1.grd -Sb${AMIN}/$ABELOW 2>/dev/null # clip minimum
        gmt grdmath modelmask.grd tmp1.grd MUL = outf_${comp}.grd 2>/dev/null # clip to TS region
        gmt grdclip tmp1.grd -Goutf_${comp}.grd -Sb${AMIN}/$ABELOW 2>/dev/null #clip minimum

        \cp outf_${comp}.grd tmp1.grd
        gmt grdmath modelmask.grd 1 SUB tmp1.grd ADD = outf_${comp}.grd 2>/dev/null
        #add grid image to ps plot
        gmt grdimage outf_${comp}.grd $att -C$BASECPT -Q -t50 -K -O >> $plot_file 2>/dev/null
        #add coastline
        gmt pscoast -A0/0/1 -N1 -N2 $att -Df -S135/205/250 -W1,black -K -O >> $plot_file
        gmt pscoast -A0/2/2 $att -Df -W1,black -K -O >> $plot_file

        #call fault plane routine
        bash ${ADDFLTPLANEDIR}/addStandardFaultPlane.sh $plot_file -R$plot_ts_region -JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} $FAULTFILE $LINE $TOPEDGE $HYPOPEN

        #add Main and subtitles
        gmt pstext $att -N -O -K -D0.0/0.35 -F+f20p,Helvetica-Bold,black+jLB+a0 << END >>  $plot_file #
$plot_x_min $plot_y_max $plot_main_title
END

        # subtitle
        gmt pstext $att -N -O -K -D0.0/0.1 -F+f+j+a0, << END >>  $plot_file
$plot_x_min $plot_y_max 14,Helvetica,black LB $plot_sub_title
$plot_x_max $plot_y_max 16,Helvetica,black RB t=$tt sec
END

        #add scale to show distance
        gmt psbasemap $att -L172.50/-43.90/${avg_ll[1]}/25.0 -Ba30mf30mWSen -K -O >> $plot_file

        # add sites
        for i in "${!plot_s_lon[@]}"; do
            add_site $i
        done

        # plot strong motion station locations
        gmt psxy "$stat_file" $att -St0.08 -G000/000/000 -W$plot_s_lin -O -K >> $plot_file
        # shift plotting origin (for 3 component plotting)
        gmt psxy -V $att -L -W5,255/255/0 -O -K -X$plot_x_shift << END >>  $plot_file 2>/dev/null
END
        #end of component
    done

    #finalize the plot (i.e. no -K)
    gmt psxy -V $att -L -W5,255/255/0 -O << END >>  $plot_file 2>/dev/null
END
    # ps -> png
    gmt ps2raster $plot_file -A -TG -E$plot_res -D$plot_png_dir

    # work on next slice
    tsfcnt=$(($tsfcnt + 1))
done

# temporary files
\rm tmp0.grd tmp1.grd modelmask.grd landmask.grd allmask.grd
\rm tmp.modelpath outf.? outf_?.grd
\rm $BASECPT
\rm gmt.conf gmt.history







wait
echo
echo Plotting Finished.
echo


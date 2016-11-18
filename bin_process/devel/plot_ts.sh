#!/usr/bin/env bash

# Created: 21 April 2016
# Purpose: Generate visualisations of timeslices (postscript and png)
# Replacing Purpose: Use bash instead of csh. Source vars from e3d.par
# Authors: Viktor Polak <viktor.polak@canterbury.ac.nz>

# USAGE:
# Execute from current directory only: "$ ./plot_ts.sh" or "$ bash plot_ts.sh"
# Optional first parameter: specify number of threads
#   Default is user interactive with autodetect, capped to 8 threads

# ISSUES:
# could validate parameters/check if folders/files exist.
# only remove .png/.ps/.bb/.x files which are actually generated/used as temp instead of *
# make executable from any location (note gmt.conf and e3d.par location)

# most variables are used from e3d.par
source e3d.par
# script is run with second parameter to indicate testing/override parameters
if [ "$2" != '' ]; then
    echo Running with test params from $2
    source "$2"
fi

# threading: use first parameter if available (prevents user interaction)
threads=$1
if [ $threads -gt 0 ] 2>/dev/null; then
    echo Using $threads threads from passed in value.
elif [ "$1" != "template" ]; then
    virtual_cores=$(lscpu | grep '^CPU(s):')
    virtual_cores=${virtual_cores##* }
    if [ $virtual_cores -gt 8 ]; then
        # performance issues with more threads (disk IO?)
        virtual_cores=8
    fi

    read -p "Run on how many threads? [$virtual_cores]: " threads
    if [ "$threads" == '' ]; then
        threads=$virtual_cores
        echo Using default, $virtual_cores threads.
    elif [ "$threads" -gt 0 ] 2>/dev/null; then
        echo Running on $threads threads.
    else
        echo Invalid input. Exiting.
        exit
    fi
fi

# threading doesn't make sense with template option
# only create template ps/raster for later use (web use)
if [ "$1" == "template" ] || [ "$1" == "template2" ]; then
    plot_purpose="$1"
    if [ ! -e "$tsbin.bin" ]; then
        echo Need to run hdf2maxgrid.py first.
        echo Cannot find $tsbin.bin.
        exit
    fi
fi

# definition of different regions to plot for
case $plot_region in
    # region to plot
    # sites to display
    # place dot at [Left Centre Right, Top Middle Bottom] of text
    # scale to show distance
    CANTERBURY)
        plot_x_min=171.75
        plot_x_max=173.00
        plot_y_min=-44.00
        plot_y_max=-43.20
        plot_sites=(Rolleston Darfield Lyttelton Akaroa Kaiapoi Rakaia Oxford)
        plot_s_pos=(RB CB LM RB LB RT LB)
        plot_s_lon=(172.3791667 172.1116667 172.7194444 172.9683333 172.6569444 172.0230556 172.1938889)
        plot_s_lat=(-43.59083333 -43.48972222 -43.60305556 -43.80361111 -43.38277778 -43.75611111 -43.29555556)
        plot_scale="-L172.50/-43.90/$modellat/25.0 -Ba30mf30mWSen"
        ;;
    WIDERCANT)
        plot_x_min=170.52
        plot_x_max=173.67
        plot_y_min=-44.4
        plot_y_max=-42.53
        plot_sites=(Rolleston Darfield Lyttelton Akaroa Kaiapoi Rakaia Oxford)
        plot_s_pos=(RB CB LM RB LB RT LB)
        plot_s_lon=(172.3791667 172.1116667 172.7194444 172.9683333 172.6569444 172.0230556 172.1938889)
        plot_s_lat=(-43.59083333 -43.48972222 -43.60305556 -43.80361111 -43.38277778 -43.75611111 -43.29555556)
        plot_scale="-L172.50/-43.90/$modellat/25.0 -Ba30mf30mWSen"
        ;;
    SOUTHISLAND)
        plot_x_min=166.0
        plot_x_max=174.5
        plot_y_min=-47.50
        plot_y_max=-40.00
        plot_sites=(Queenstown Dunedin Tekapo Timaru Christchurch Haast Greymouth Westport Kaikoura Nelson Blenheim)
        plot_s_pos=(LM LM LM LM LM RM RM RM CT CB CT)
        plot_s_lon=(168.6680556 170.3794444 170.4794444 171.2430556 172.6347222 169.0405556 171.2063889 171.5997222 173.6802778 173.2838889 173.9569444)
        plot_s_lat=(-45.0300000 -45.8644444 -44.0069444 -44.3958333 -43.5313888 -43.8808333 -42.4502777 -41.7575000 -42.4038888 -41.2761111 -41.5138888)
        plot_scale="-L173/-47/$modellat/50.0 -Ba60mf60mWSen"
        ;;
    MIDNZ)
        plot_x_min=168.2
        plot_x_max=177.9
        plot_y_min=-45.7
        plot_y_max=-37.85
        plot_sites=(Queenstown Tekapo Timaru Christchurch Haast Greymouth Westport Kaikoura Nelson Blenheim Wellington Palmerston\ North Masterton Napier New\ Plymouth Taupo Rotorua)
        plot_s_pos=(LM LM LM LM RM RM RM LM CB LM RM RM LM LM RM LM LM)
        plot_s_lon=(168.6680556 170.4794444 171.2430556 172.6347222 169.0405556 171.2063889 171.5997222 173.6802778 173.2838889 173.9569444 174.777222 175.611667 175.658333 176.916667 174.083333 176.069400 176.251389)
        plot_s_lat=(-45.0300000 -44.0069444 -44.3958333 -43.5313888 -43.8808333 -42.4502777 -41.7575000 -42.4038888 -41.2761111 -41.5138888 -41.288889 -40.355000 -40.952778 -39.483333 -39.066667 -38.6875 -38.137778)
        plot_scale="-L176/-45/$modellat/100.0 -Ba60mf60mWSen"
        ;;
    *)
        echo Plotting Region Not Understood
        exit
        ;;
esac
# avg lon/lat (midpoint of plotting region)
avg_ll=(`echo $plot_x_min $plot_x_max $plot_y_min $plot_y_max | \
        gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`)
# GMT format boundaries
plot_ll_region=$plot_x_min/$plot_x_max/$plot_y_min/$plot_y_max
# calculate distance on figure (page size) for plotting region
# only image accurate if using a straight lat/lon projection
# inch / degree
ipd=0.6
x_size=( $(echo $plot_x_min $plot_x_max $ipd | \
        gawk '{len = ($2 - $1) * $3; print 0, len / 2, len}') )
y_size=( $(echo $plot_y_min $plot_y_max $ipd | \
        gawk '{len = ($2 - $1) * $3; print 0, len / 2, len}') )

# output dirs
mkdir -p $plot_ps_dir $plot_png_dir
rm "$plot_ps_dir"/* "$plot_png_dir"/* 2>/dev/null
# temporary working directories for gmt are within here
# gmt process create temp files which clash between process
gmt_temp="$sim_dir/gmt_wd"
mkdir -p "$gmt_temp"

# make temp file with simulation boundaries
if [ -e "sim.modelpath" ]; then
    \rm sim.modelpath
fi
# have to close box (copy first corner to end)
for corner in c1 c2 c3 c4 c1; do
    grep "$corner= " $vel_mod_params | gawk ' { print $2, $3 } ' >> sim.modelpath
    # only plot within simulation domain for adding to google maps
    if [ "$plot_purpose" == "template" ]; then
        case "$corner" in
            c1)
                west=$(grep "$corner= " $vel_mod_params | gawk ' { print $2 } ')
                ;;
            c2)
                north=$(grep "$corner= " $vel_mod_params | gawk ' { print $3 } ')
                ;;
            c3)
                east=$(grep "$corner= " $vel_mod_params | gawk ' { print $2 } ')
                ;;
            c4)
                south=$(grep "$corner= " $vel_mod_params | gawk ' { print $3 } ')
                ;;
        esac
    fi
done

# create color palette for plotting the topography
base_cpt=y2r_brown.cpt

add_sites() {
    psfile=$1

    for i in "${!plot_s_lon[@]}"; do
        # location as a circle. black outline, semi-transparent grey fill.
        echo ${plot_s_lon[$i]} ${plot_s_lat[$i]} | \
                psxy -R -J -Sc0.10 -G220/220/220@50 -W0.8,black -O -K >> "$1"

        # location name
        echo ${plot_s_lon[$i]} ${plot_s_lat[$i]} ${plot_s_pos[$i]} ${plot_sites[$i]} | \
                pstext -R -J -N -O -K -Dj0.08/0.08 -F+j+f10,Helvetica,black+a0 >>  "$1"
    done
}

add_source() {
    # add the fault plane or beachball
    if [ "$plot_type" == "finitefault" ]; then
        # plot fault plane
        bash ${plot_fault_add_plane} "$1" \
            -R$plot_ll_region -J \
            $fault_file $plot_fault_line $plot_fault_top_edge $plot_fault_hyp_open
    else
        # plot beach ball (variable must be surrounded by quotes in sourced file)
        echo $plot_beachball | \
                gmt psmeca -P -J -R -Sc0.15 -Ggreen -O -K >> "$1"
    fi
}

finalise_png() {
    # finalises PS and creates PNG
    psfile=$1
    outdir=$2
    outres=$3
    # finalize postscript (i.e. no -K)
    psxy -J -R -O -T >>  "$psfile" 2>/dev/null
    # ps -> png
    # ps2raster deprecated in 5.2, replaced by psconvert
    if [ -x "$(which psconvert 2>/dev/null)" ]; then
        psconvert "$psfile" -A -TG -E"$outres" -D"$outdir"
    else
        ps2raster "$psfile" -A -TG -E"$outres" -D"$outdir"
    fi
}

clean_temp_files() {
    # temporary files (global)
    \rm modelmask.grd
    \rm $base_cpt
    \rm gmt.conf gmt.history

    # gmt process working directories
    rm -r "$gmt_temp"
}

sig_int_received() {
    # restore INT behaviour
    trap "" INT
    echo
    echo
    echo Caught interrupt. Waiting for processes and cleaning directory...
    echo
    echo
    # kill children
    pkill -P $$

    clean_temp_files
    echo Process exiting.
    exit
}
# enable killing of all subprocesses/cleaning temp files on CTRL-C (interrupt)
trap "sig_int_received" INT

# color palette for velocity TODO: make into parameter
# https://www.soest.hawaii.edu/gmt/gmt/html/images/GMT_RGBchart_a4.png
#if [ "$absmax" -eq 1 ] || [ "$plot_purpose" == "template" ]; then
    cpt=hot
    min=0
    extra='-I'
    scale='ground velocity (cm/s)'
    seismoline=blue
#else
#    cpt=polar
#    min=-$plot_topo_a_max
#    extra=''
#    scale='ground velocity east (cm/s)'
#    seismoline=darkgreen
#fi
makecpt -C$cpt $extra -T$min/$plot_topo_a_max/$plot_topo_a_inc -A50 > $base_cpt

# set all plotting defaults to use
gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i PS_MEDIA = A2

# create mask for simulation domain
grdmask sim.modelpath_hr -NNaN/0/1 -R$plot_ll_region -I$plot_dx/$plot_dy -Gmodelmask.grd

###################### BEGIN TEMPLATE ##########################
echo Creating PS Template...
plot_file_template=$gmt_temp/plot_template.ps
# draw background given main image offsets (inches)
# left and bottom are exact when using uniform longitude projection
# top is fiddly as height of image depends on projection distortion
# additional white background is cropped (PNG), so is anything below/left of it
left=1.5
bottom=1.4
right=0.8
top=0.8
# counter shift ($1) needs to be multiplied by 2 in -JX parameter
# the distance is wrong after *2 but doesn't matter with -JX crop
length=$(echo $left $right ${x_size[2]} | gawk '{print $1 * 2 + $2 + $3}')
height=$(echo $bottom $top ${y_size[2]} | gawk '{print $1 * 2 + $2 + $3 * 1.384}')
# draw background after shifting plotting origin
echo -e "$length 0\n$length $height\n0 $height\n0 0" | \
        psxy -G255/255/255 -JX$length/$height -R0/$length/0/$height \
                -K -Xa-$left -Ya-$bottom > "$plot_file_template"
# switch origin back and set projection/region
psxy -T -Jm${avg_ll[0]}/$ipd -R$plot_ll_region -X$left -Y$bottom -O -K >> "$plot_file_template"


# add land (usually fully covered by topography layer)
pscoast -R$plot_ll_region -J -Ba60mf30mWSen -Df -G200/200/200 -K -O >> "$plot_file_template"
# topography palette re-scaling
makecpt -Cgray -T-5000/3000/1 > land.cpt
# make topography given topo and illumination file
grdimage $plot_topo_file $plot_topo_illu -Cland.cpt -Q -R -J -K -O >> "$plot_file_template"
rm land.cpt
# add coastline outline on top of topography for full line
#pscoast -R -J -Df -W1p,black -K -O >> "$plot_file_template"
# add oceans/coastline
pscoast -A0/0/2 -Na/1.0p,black -R -J -Df -S95/165/250 -K -O >> "$plot_file_template"
# add lakes/coastline
pscoast -A0/2/2 -R -J -Df -S95/165/250 -K -O >> "$plot_file_template"
# simulation domain
gmt psxy sim.modelpath -R -J -W0.4p,black,- -L -O -K >> "$plot_file_template"

# set the color scale
# -Efb (forground/background triangles) -Dcentre/ydisplacement/length/height -Ba<major tick>f<minor tick>
psscale -C$base_cpt -Ef -D${x_size[1]}/-0.5/4.0/0.15h -K -O -Ba${plot_topo_a_inc}f${plot_topo_a_inc}:"ground motion (cm/s)": >> "$plot_file_template"
# main title
echo ${avg_ll[0]} $plot_y_max $plot_title | \
        pstext -R -J -N -O -K -D0/0.4 \
                -F+f20p,Helvetica,black+jCB >>  "$plot_file_template"
# subtitle
echo $plot_x_min $plot_y_max 16,Helvetica,black LB "$plot_sub_title" | \
        pstext -R -J -N -O -K -D0.0/0.1 -F+f+j+a0, >>  "$plot_file_template"

if [ "$plot_purpose" == "template" ]; then
    # subtitle
    echo $plot_x_max $plot_y_max 16,Helvetica,black RB t=180.0 sec | \
            pstext -R -J -N -O -K -D0.0/0.1 -F+f+j+a0, >>  "$plot_file_template"

    # include maxgrid on map
    # create ground motion intensity surface from MAXGRID
    if [ "$swap_bytes" -eq 1 ]; then
        xyz2grd $tsbin.bin -S$tsbin.native.bin -V -Zf 2>/dev/null
    else
        cp $tsbin.bin $tsbin.native.bin
    fi
    surface $tsbin.native.bin -Gmax.grd -I$plot_dx/$plot_dy \
            -R$plot_ll_region -T0.0 -bi3f # 2>/dev/null
    rm $tsbin.native.bin
    # crop to simulation domain (multiply by mask of 0 or 1, outside = 0)
    grdmath max.grd modelmask.grd MUL = max.grd 2>/dev/null
    # clip minimum (values below cutoff = NaN, not displayed, clear)
    grdclip max.grd -Gmax.grd \
            -Sb${plot_topo_a_min}/NaN #2>/dev/null

    # clippath for land
    #pscoast -R -J -Df -Gc -K -O >> "$plot_file_template"
    # add resulting overlay image to plot
    grdimage max.grd -t40 -R -J -C$base_cpt -Q -K -O >> "$plot_file_template" 2>/dev/null
    # clear clippath (crop ocean area)
    #pscoast -R -J -O -K -Q >> "$plot_file_template"

    # add sites
    add_sites "$plot_file_template"
    # finite fault or beachball
    add_source "$plot_file_template"
    rm sim.modelpath
    # add check points (used when aligning on overlay)
    #    psxy -R -J -Sc0.005 -G000/255/255 -W1.0 -O -K << merda >> "$plot_file_template"
    #171.4688139 -41.7448556
    #169.876852777777778 -42.424077777777778
    #168.488294444444444445 -46.515325
    #172.6899 -43.78415555555555555556
    #171.7632361111111 -43.25080555555555
    #merda
    mkdir -p WEB
    cp "$plot_file_template" WEB/gmt-template.ps
    cp gmt.conf WEB/
    cp gmt.history WEB/
    finalise_png "$plot_file_template" "WEB" "600"
    clean_temp_files
    # create kmz
    cd WEB
    echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > "doc.kml"
    echo "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">" >> "doc.kml"
    # color is alpha, r,g,b
    echo -e "<GroundOverlay>\n\t<name>Max GM Overlay</name>\n\t<color>80ffffff</color>\n\t<Icon>\n\t\t<href>plot_template.png</href>\n\t\t<viewBoundScale>0.75</viewBoundScale>\n\t</Icon>\n\t<altitude>22206</altitude>\n\t<LatLonBox>\n\t\t<north>$north</north>\n\t\t<south>$south</south>\n\t\t<east>$east</east>\n\t\t<west>$west</west>\n\t</LatLonBox>\n</GroundOverlay>\n</kml>" >> "doc.kml"
    zip motion.kmz plot_template.png doc.kml
    exit
fi
echo Template Complete.
####################### END TEMPLATE ###########################

render_slice() {
    # $1 is $tsfcnt
    # $2 is $tsspot
    # $3 is $tsfnum

    # used to display t=tt sec in subtitle
    tt=`echo $2 $plot_orig_dt | gawk '{printf "%.2f\n",$1*$2;}'`

    # outputs
    plot_file=$plot_ps_dir/ts-str${3}.ps
    png_file=$plot_png_dir/ts-str${3}.png
    cp "$plot_file_template" "$plot_file"

    ### TEMP moved before ground motion
    # plot strong motion station locations
    psxy "$stat_file" -R -J -St0.08 -W0.5p,black -O -K >> "$plot_file"
    # add coastline outline on top of topography for full line
    pscoast -R -J -Df -W0.3,black -K -O >> "$plot_file"

    # basename for timeslice inputs
    outf=`echo $ts_out_prefix $1 | gawk '{printf "%s_ts%.4d\n",$1,$2;}'`

    if [ "$swap_bytes" -eq 1 ]; then
        # all components are always separate by new get_ts version
        # if for some reason using old data:
        #    set absmax=0, comment last 2 lines in both blocks, +more
        xyz2grd ${outf}.0 -Soutf_${3}.0 -V -Zf 2>/dev/null
        #xyz2grd ${outf}.1 -Soutf_${3}.1 -V -Zf 2>/dev/null
        #xyz2grd ${outf}.2 -Soutf_${3}.2 -V -Zf 2>/dev/null
    else
        \cp ${outf}.0 outf_${3}.0
        #\cp ${outf}.1 outf_${3}.1
        #\cp ${outf}.2 outf_${3}.2
    fi

    # create ground motion intensity surface from the TSlice output
    # -bi for binary input of 3 columns of floats
    surface outf_${3}.0 -Gtmp_${3}.grd -I$plot_dx/$plot_dy \
            -R -T0.0 -bi3f 2>/dev/null
    if [ "$absmax" -eq 1 ]; then
        # velocity = SQRT(X^2 + Y^2 + Z^2)
        surface outf_${3}.1 -Gtmp_${3}.1.grd -I$plot_dx/$plot_dy \
                -R$plot_ll_region -T0.0 -bi3f 2>/dev/null
        surface outf_${3}.2 -Gtmp_${3}.2.grd -I$plot_dx/$plot_dy \
                -R$plot_ll_region -T0.0 -bi3f 2>/dev/null
        grdmath tmp_${3}.grd SQR tmp_${3}.1.grd SQR ADD tmp_${3}.2.grd SQR ADD SQRT = tmp_${3}.grd 2>/dev/null
        rm tmp_${3}.1.grd tmp_${3}.2.grd outf_${3}.1 outf_${3}.2
    fi
    # crop to simulation domain (multiply by mask of 0 or 1, outside = 0)
    grdmath tmp_${3}.grd modelmask.grd MUL = tmp_${3}.grd 2>/dev/null
    # clip minimum (values below cutoff = NaN, not displayed, clear)
    grdclip tmp_${3}.grd -Gtmp_${3}_P.grd \
            -Sb${plot_topo_a_min}/NaN 2>/dev/null
    #if [ "$absmax" -eq 0 ]; then
    #    # also cutoff for the lower section
    #    grdclip tmp_${3}.grd -Gtmp_${3}_N.grd \
    #            -Sa-${plot_topo_a_min}/NaN 2>/dev/null
    #    # clip max, colour scale only goes down to > -$plot_topo_a_max
    #    # if the scale is inverted, have to clip min instead
    #    grdmath $plot_topo_a_max 1 SUB tmp_${3}_P.grd MIN tmp_${3}_N.grd AND = tmp_${3}_P.grd 2>/dev/null
    #    rm tmp_${3}_N.grd
    #fi

    # clippath for land
    #pscoast -R -J -Df -Gc -K -O >> "$plot_file"
    # add resulting overlay image to plot
    grdimage tmp_${3}_P.grd -R -J -C$base_cpt -Q -t40 -K -O >> "$plot_file" 2>/dev/null
    # clear clippath (crop ocean area)
    #pscoast -R -J -O -K -Q >> "$plot_file"

    # remove temporary input (potentially byte swapped) and grid file
    rm outf_${3}.0 tmp_${3}.grd tmp_${3}_P.grd

    # ADDFAULTPLANE.SH MAKES TEMP FILES WHICH INTERFERE (SAME NAME)
    # CD INTO TEMP DIR (REQUIRES ABS PATHS)
    # TODO: fix addfaultplane or implement it here
    gmt_proc_wd=$(mktemp -d -p "$gmt_temp" GMT.XXXXXXXX)
    cd "$gmt_proc_wd"
    #cp "$sim_dir/gmt.conf" ./
    #cp "$sim_dir/gmt.history" ./
    # for testing, could use relative paths instead
    cp ../../gmt.conf ./
    cp ../../gmt.history ./
    plot_file=../../$plot_file
    plot_png_dir=../../Png

    # subtitle part 2 (dynamic)
    # must precede psmeca as font configuration history required
    echo $plot_x_max $plot_y_max 16,Helvetica,black RB t=$tt sec | \
            pstext -R -J -N -O -K -D0.0/0.1 -F+f+j+a0, >>  "$plot_file"

    # finite fault or beachball
    add_source "$plot_file"

    # scale to show distance
    psbasemap -R -J $plot_scale -K -O >> "$plot_file"

    # add sites
    add_sites "$plot_file"

    # add seismograms
    # --no-group-separator only for GNU grep, if other grep, use -v '^--$'
    if [ -e "../../gmt-seismo.xy" ]; then
        if [ "$plot_seismo_style" == "SCEC" ]; then
            pattern="^>TS$2\s"
        else
            pattern='^>'
        fi
        psxy -N -R -J -W1.5p,$seismoline -O -K << END >> "$plot_file"
$(cat ../../gmt-seismo.xy | grep -e "$pattern" -A $(($2 + 1)) --no-group-separator)
END
    fi

    # finalises PostScript and converts to PNG
    finalise_png "$plot_file" "$plot_png_dir" "$plot_res"

    # remove own working directory
    cd - >/dev/null
    rm -r "$gmt_proc_wd"
}

# start loop for all time slices
tsfcnt=$ts_start # variable used in process thread
while [ $tsfcnt -lt $ts_total ]; do
    # current TS index - used in process thread
    tsspot=$(($tsfcnt * $ts_inc))
    # output number (4 decimal places) - used in process thread
    tsfnum=`echo $tsfcnt | gawk '{printf "%.4d\n",$1;}'`

    # wait until running threads are less than max
    while [ $(jobs | wc -l ) -ge $threads ]; do
        sleep 0.5
        jobs >/dev/null
    done

    # call the worker process
    render_slice $tsfcnt $tsspot $tsfnum &
    pid=$!
    echo Rendering timeslice $(($tsfcnt + 1))/$ts_total SEQ: $tsfcnt TS: $tsspot PID: $pid

    # work on next slice
    tsfcnt=$(($tsfcnt + 1))
done

wait
clean_temp_files
# restore INT behaviour
trap "" INT

echo
echo Plotting Finished.
echo


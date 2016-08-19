#!/usr/bin/env bash

# Created: 21 April 2016
# Purpose: Generate visualisations of timeslices (postscript and png)
# Replacing: plot_ts_bluefern.csh
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

# COMPATIBILITY/PORTABILITY:
# bash compatible, should work with ksh too.
# Bashisms: arrays, '==' in tests, read (internal command) parameters
# Linuxisms: 'lscpu' output, only tested with GNU AWK

# following variables are used from e3d.par:
#   plot_main_title, plot_sub_title, plot_option, wcc_prog_dir,
#   stat_file, vel_mod_params_dir, h, dt_ts, dx_ts, dy_ts, dz_ts, ts_start, ts_inc,
#   swap_bytes, lonlat_out, scale, plot_sites, plot_s_pos, plot_s_lon, plot_s_lat,
#   plot_s_sym, plot_s_fil, plot_s_lin, plot_x_org, plot_y_org, plot_x_inch, plot_x_shift,
#   plot_region, plot_dx, plot_dy, plot_palette,
#   ts_file, ts_out_prefix, plot_ps_dir, plot_png_dir, plot_res, plot_orig_dt, plot_comps,
#   global_root, sim_dir, plot_topo_file, plot_topo_illu, plot_topo_a_min, plot_topo_a_inc,
#   plot_topo_a_max, plot_topo_a_below, plot_fault_{add_plane,line,top_edge,hyp_open}, fault_file,
#   modellat
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
else # invalid or no input
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

# definition of different regions to plot for
case $plot_region in
    # region to plot
    # sites to display
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
    SOUTHISLAND)
        plot_x_min=166.0
        plot_x_max=174.5
        plot_y_min=-47.50
        plot_y_max=-40.00
        plot_sites=(Queenstown Dunedin Tekapo Timaru Christchurch Haast Greymouth Westport Kaikoura Nelson Blenheim)
        plot_s_pos=(LT LM LM LM LM RM RM RM LM CB LM)
        plot_s_lon=(168.6680556 170.3794444 170.4794444 171.2430556 172.6347222 169.0405556 171.2063889 171.5997222 173.6802778 173.2838889 173.9569444)
        plot_s_lat=(-45.0300000 -45.8644444 -44.0069444 -44.3958333 -43.5313888 -43.8808333 -42.4502777 -41.7575000 -42.4038888 -41.2761111 -41.5138888)
        plot_scale="-L173/-47/$modellat/50.0 -Ba60mf60mWSen"
        ;;
    *)
        echo Plotting Region Not Understood
        exit
        ;;
esac
# common region related properties
plot_ll_region=$plot_x_min/$plot_x_max/$plot_y_min/$plot_y_max
plot_s_sym="c0.10"
plot_s_fil="220/220/220"
plot_s_lin="1,000/000/000"


# used by 'ts2xyz' bin, not used with plot_option=2
grid_file="${vel_mod_params_dir}/gridout_nz01-h${h}"
model_params="${vel_mod_params_dir}/model_params_nz01-h${h}" # to get location corners

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
    grep "$corner= " $model_params | gawk ' { print $2, $3 } ' >> sim.modelpath
done
# create mask for simulation domain
grdmask sim.modelpath -R$plot_ll_region -I$plot_dx/$plot_dy -Gmodelmask.grd
# sim.modelpath still used to plot simulation domain

# create color palette for plotting the topography
base_cpt=y2r_brown.cpt

# avg Lon/Lat (midpoint of the TSlice image for the geo projection)
avg_ll=(`echo $plot_x_min $plot_x_max $plot_y_min $plot_y_max | gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`)
# plot projection and region in GMT format for ease of use later
att="-JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} -R${plot_ll_region}"

add_site() {
    # $1 is the site index
    # $2 is the plot file

    # plot location as a point
    psxy $att -S$plot_s_sym -G$plot_s_fil -W$plot_s_lin -O -K << END >> "$2"
${plot_s_lon[$1]} ${plot_s_lat[$1]}
END

    # add location name
    pstext $att -N -O -K -Dj0.05/0.05 -F+j+f12,Helvetica,black+a0 << END >>  "$2"
${plot_s_lon[$1]} ${plot_s_lat[$1]} ${plot_s_pos[$1]} ${plot_sites[$1]}
END
}

finalise_png() {
    # finalize postscript (i.e. no -K)
    psxy -V $att -L -W5,255/255/0 -O << END >>  "$1" 2>/dev/null
END
    # ps -> png
    # ps2raster deprecated in 5.2, replaced by psconvert
    if [ -x "$(which psconvert 2>/dev/null)" ]; then
        psconvert "$1" -A -TG -E$plot_res -D$plot_png_dir
    else
        ps2raster "$1" -A -TG -E$plot_res -D$plot_png_dir
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

# color palette for velocity
#makecpt -Chot -I -T0/$plot_topo_a_max/$plot_topo_a_inc -A50 > $base_cpt
makecpt -Chot -I -T0/$plot_topo_a_max/0.01 -A50 > $base_cpt

# set all plotting defaults to use
gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i

###################### BEGIN TEMPLATE ##########################
echo Creating PS Template...
plot_file_template=$gmt_temp/plot_template.ps
# specify plot and panel size (defaults 8.5 x 11)
edge_colour=255/255/255 #180/180/180 = grey ; 255/255/255=white
psxy -JX8.5/21 -R0/8.5/0/11 -L -G${edge_colour} -X0 -Y0 -K << END > "$plot_file_template" #-W0/180/180/180
0.3 1.0
0.3 7.8
6.5 7.8
6.5 1.0
END
# set the color scale
psscale -C$base_cpt -Ef -D3.0/2.0/2.5/0.15h -K -O -Ba${plot_topo_a_inc}f${plot_topo_a_inc}:"ground velocity (cm/s)": >> "$plot_file_template"
# specify the X and Y offsets` for plotting (I dont really understand this yet)
psxy -V $att -L  -K -O -X$plot_x_org -Y$plot_y_org << END >> "$plot_file_template" 2>/dev/null #-W5/255/255/0
END
# try a different version of plotting
# clippath for land
pscoast $att -Df -Gc -K -O >> "$plot_file_template"
# land
grdimage $plot_topo_file $plot_topo_illu $plot_palette $att -K -O >> "$plot_file_template"
# clear clippath
pscoast -R -J -O -K -Q >> "$plot_file_template"
# add urban areas
URBANDIR=${global_root}/PlottingData/sourcesAndStrongMotionStations
psxy ${URBANDIR}/ChchUrbanBoundary.xy $att -G160/160/160 -W0.5p -O -K >> "$plot_file_template"
# add oceans/coastline
pscoast -A0/0/1 -N1 -N2 $att -Df -S135/205/250 -W1,black -K -O >> "$plot_file_template"
# add lakes/coastline
pscoast -A0/2/2 $att -Df -S135/205/250 -W1,black -K -O >> "$plot_file_template"
# simulation domain
gmt psxy sim.modelpath $att -W1.5p,black,- -L -O -K >> "$plot_file_template"
rm sim.modelpath
# main title
pstext $att -N -O -K -D0.0/0.35 \
        -F+f20p,Helvetica-Bold,black+jLB+a0 << END >>  "$plot_file_template"
$plot_x_min $plot_y_max $plot_main_title
END
# subtitle part 1 (static)
pstext $att -N -O -K -D0.0/0.1 -F+f+j+a0, << END >>  "$plot_file_template"
$plot_x_min $plot_y_max 14,Helvetica,black LB $plot_sub_title
END
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

    if [ "$plot_option" -eq 1 ]; then
        # use 'ts2xyz.exe. to get TSlice outout in xyz format
        # XXX: WARNING: outfile has been changed from 'outf' to 'outf_${3}'
        #   not sure what is happening with it (appending or needing to be fresh)
        ${wcc_prog_dir}/ts2xyz infile=$ts_file outfile=outf_${3} swap_bytes=$swap_bytes \
                gridfile=$grid_file xyts=1 scale=$scale ts=$2 trv=0 \
                dxts=$dxts dyts=$dyts dzts=1 absmax=1 \
                read_header=1 outbin=1 lonlat=$lonlat_out geoproj=1
    elif [ "$plot_option" -eq 2 ]; then
        outf=`echo $ts_out_prefix $1 | gawk '{printf "%s_ts%.4d\n",$1,$2;}'`
    fi

    # currently components are not looped over, as using ABSMAX=1
    # loop over the different components
    for comp in "${plot_comps[@]}"; do
        if [ "$plot_option" -eq 2 ]; then
            if [ "$swap_bytes" -eq 1 ]; then
                # get file in correct format - BB added
                xyz2grd ${outf}.${comp} -Soutf_${3}.${comp} -V -Zf
            elif [ "$swap_bytes" -eq 0 ]; then
                \cp ${outf}.${comp} outf_${3}.${comp}
            fi
        fi

        # create ground motion intensity surface from the TSlice output
        # -bi for binary input of 3 columns of floats
        surface outf_${3}.${comp} -Gtmp_${3}.grd -I$plot_dx/$plot_dy \
                -R$plot_ll_region -T0.0 -bi3f 2>/dev/null
        # crop to simulation domain (multiply by mask of 0 or 1, outside = 0)
        grdmath tmp_${3}.grd modelmask.grd MUL = tmp_${3}.grd 2>/dev/null
        # clip minimum (values below cutoff = NaN, not displayed, clear)
        grdclip tmp_${3}.grd -Gtmp_${3}.grd \
                -Sb${plot_topo_a_min}/NaN 2>/dev/null
        # add resulting overlay image to plot
        grdimage tmp_${3}.grd $att -C$base_cpt -Q -t50 -K -O >> "$plot_file" 2>/dev/null
        # remove temporary input (potentially byte swapped) and grid file
        rm outf_${3}.${comp} tmp_${3}.grd

        #psxy  Roads.gmt >> "$plot_file"
        # add resulting overlay image to plot
        #grdimage roads.grd $att -Q -K -O >> "$plot_file"

        # ADDFAULTPLANE.SH MAKES TEMP FILES WHICH INTERFERE (SAME NAME)
        # CD INTO TEMP DIR (REQUIRES ABS PATHS)
        # TODO: fix addfaultplane or implement it here
        gmt_proc_wd=$(mktemp -d -p "$gmt_temp" GMT.XXXXXXXX)
        cd "$gmt_proc_wd"
        #cp "$sim_dir/gmt.conf" ./
        # for testing, could use relative path instead
        cp ../../gmt.conf ./

        # subtitle part 2 (dynamic)
        # must precede psmeca as font configuration history required?
        pstext $att -N -O -K -D0.0/0.1 -F+f+j+a0, << END >>  "$plot_file"
$plot_x_max $plot_y_max 16,Helvetica,black RB t=$tt sec
END

        # add the fault plane or beachball
        if [ "$plot_type" == "finitefault" ]; then
            # plot fault plane
            bash ${plot_fault_add_plane} "$plot_file" \
                -R$plot_ll_region -JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} \
                $fault_file $plot_fault_line $plot_fault_top_edge $plot_fault_hyp_open
        else
            # plot beach ball (variable must be surrounded by quotes in sourced file)
            gmt psmeca -P -J -R -Sc0.15 -Ggreen -O -K << EOF >> "$plot_file"
$plot_beachball
EOF
        fi

        # scale to show distance
        psbasemap $att $plot_scale -K -O >> "$plot_file"

        # add sites
        for i in "${!plot_s_lon[@]}"; do
            add_site "$i" "$plot_file"
        done

        # plot strong motion station locations
        psxy "$stat_file" $att -St0.08 -G000/000/000 -W$plot_s_lin -O -K >> "$plot_file"

        # add seismograms
        psxy ../../station.xy $att -W1.5p,red -O -K >> "$plot_file"

        # shift plotting origin (for 3 component plotting)
        psxy -V $att -L -W5,255/255/0 -O -K -X$plot_x_shift << END >>  "$plot_file" 2>/dev/null
END
    done

    # finalises PostScript and converts to PNG
    finalise_png "$plot_file"

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


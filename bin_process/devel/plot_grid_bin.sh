#!/usr/bin/env bash

# Created: 10 May 2016
# Purpose: Generate visualisation of max ground motion (postscript and png)
# Authors: Viktor Polak <viktor.polak@canterbury.ac.nz>

# USAGE:

# ISSUES:

# following variables are used from e3d.par:
#   plot_main_title, plot_sub_title, plot_option, wcc_prog_dir,
#   stat_file, vel_mod_params_dir, h, dt_ts, dx_ts, dy_ts, dz_ts, ts_start, ts_inc,
#   swap_bytes, lonlat_out, scale, plot_sites, plot_s_pos, plot_s_lon, plot_s_lat,
#   plot_s_sym, plot_s_fil, plot_s_lin, plot_x_org, plot_y_org, plot_x_inch, plot_x_shift,
#   plot_x_min, plot_x_max, plot_y_min, plot_y_max, plot_region, plot_ts_x_min, plot_ts_x_max,
#   plot_ts_y_min, plot_ts_y_max, plot_region, plot_ts_region, plot_dx, plot_dy, plot_palette,
#   ts_file, ts_out_prefix, plot_ps_dir, plot_png_dir, plot_res, plot_orig_dt, plot_comps,
#   global_root, sim_dir, plot_topo_file, plot_topo_illu, plot_topo_a_min, plot_topo_a_inc,
#   plot_topo_a_max, plot_topo_a_below, plot_fault_{add_plane,line,top_edge,hyp_open}, fault_file,
#   
source e3d.par

# outputs
plot_file=motion-map.ps
png_file=motion-map.png

# output containing lon, lat, max_val in 32bit float binary
# TODO: add to e3d.par or make python version
input='stations_max.bin'

# used by 'ts2xyz' bin, not used with plot_option=2
grid_file="${vel_mod_params_dir}/gridout_nz01-h${h}"
model_params="${vel_mod_params_dir}/model_params_nz01-h${h}" # to get location corners

# make temp file with the model coord boundaries
if [ -e "tmp.modelpath" ]; then
    \rm tmp.modelpath
fi
for corner in c1 c2 c3 c4; do
    grep "$corner= " $model_params | gawk ' { print $2, $3 } ' >> tmp.modelpath
done

# create masks for later plotting of different layers
gmt grdlandmask -R$plot_ts_region -Dh -I$plot_dx/$plot_dy -Glandmask.grd
gmt grdmask tmp.modelpath -R$plot_ts_region -I$plot_dx/$plot_dy -Gmodelmask.grd
gmt grdmath landmask.grd modelmask.grd MUL = allmask.grd
\rm tmp.modelpath

# create color palette for plotting the topography
base_cpt=y2r_hot.cpt

# avg Lon/Lat (midpoint of the TSlice image for the geo projection)
avg_ll=(`echo $plot_x_min $plot_x_max $plot_y_min $plot_y_max | gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`)
# plot projection and region in GMT format for ease of use later
att="-JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} -R${plot_region}"

# color palette for velocity, execute gmt makecpt for more info
# seis: R-O-Y-G-B -I(nverted) -T(min)/(max)/(inc)
gmt makecpt -Chot -I -T0/80/1 -A50 > $base_cpt

# set all plotting defaults to use
gmt gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i

# specify plot and panel size (defaults 8.5 x 11)
edge_colour=255/255/255 #180/180/180 = grey ; 255/255/255=white
gmt psxy -JX8.5/11 -R0/8.5/0/11 -L -G${edge_colour} -X0 -Y0 -K << END > "$plot_file" #-W0/180/180/180
0.3 1.0
0.3 7.8
6.5 7.8
6.5 1.0
END
# set the color scale
gmt psscale -C$base_cpt -Ef -D3.0/2.0/2.5/0.15h -K -O -Ba10f$10:"ground velocity (cm/s)": >> "$plot_file"
# specify the X and Y offsets for plotting (I dont really understand this yet)
gmt psxy -V $att -L  -K -O -X$plot_x_org -Y$plot_y_org << END >> "$plot_file" 2>/dev/null #-W5/255/255/0
END
# try a different version of plotting
# clippath for land
gmt pscoast $att -Df -Gc -K -O >> "$plot_file"
# land
gmt grdimage $plot_topo_file $plot_topo_illu $plot_palette $att -K -O >> "$plot_file"
# clear clippath
gmt pscoast -R -J -O -K -Q >> "$plot_file"
# add urban areas
URBANDIR=${global_root}/PlottingData/sourcesAndStrongMotionStations
gmt psxy ${URBANDIR}/ChchUrbanBoundary.xy $att -G160/160/160 -W0.5p -O -K >> "$plot_file"

# create ground motion intensity surface
gmt surface $input -Gmotion.surface -I$plot_dx/$plot_dy \
        -R$plot_ts_region -T0.0 -bi3f
# clip minimum (prevents scattered spots)
gmt grdclip motion.surface -Gmotion.surface \
        -Sb${plot_topo_a_min}/$plot_topo_a_below
# add surface as image to plot
gmt grdimage motion.surface $att -C$base_cpt -Q -t50 -K -O >> "$plot_file"

# add coastline
gmt pscoast -A0/0/1 -N1 -N2 $att -Df -S135/205/250 -W1,black -K -O >> "$plot_file"
gmt pscoast -A0/2/2 $att -Df -W1,black -K -O >> "$plot_file"

# main title
gmt pstext $att -N -O -K -D0.0/0.35 \
        -F+f20p,Helvetica-Bold,black+jLB+a0 << END >>  "$plot_file" #
$plot_x_min $plot_y_max $plot_main_title
END
# subtitle
gmt pstext $att -N -O -K -D0.0/0.1 -F+f+j+a0, << END >>  "$plot_file"
$plot_x_min $plot_y_max 14,Helvetica,black LB $plot_sub_title
$plot_x_max $plot_y_max 16,Helvetica,black RB peak velocity
END

# fault plane routine
bash ${plot_fault_add_plane} "$plot_file" \
        -R$plot_ts_region -JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} \
        $fault_file $plot_fault_line $plot_fault_top_edge $plot_fault_hyp_open

# scale to show distance
gmt psbasemap $att -L172.50/-43.90/${avg_ll[1]}/25.0 -Ba30mf30mWSen -K -O >> "$plot_file"

# finalize postscript (i.e. no -K)
gmt psxy -V $att -L -W5,255/255/0 -O << END >>  "$plot_file" 2>/dev/null
END
# ps -> png, takes a while with high resolution, 4800 works, 9600 is too high
echo Postscript complete, creating PNG...
gmt ps2raster "$plot_file" -A -TG -E600 -D./


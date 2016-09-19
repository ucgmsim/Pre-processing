#!/usr/bin/env bash
# for adding maximum motion grid overlay on google maps

source e3d.par

if [ ! -e "statgrid_max.bin" ]; then
    echo Need to run hdf2maxgrid.py first.
    echo Cannot find statgrid_max.bin.
    exit
fi

# output dirs
mkdir -p $plot_ps_dir $plot_png_dir
rm "$plot_ps_dir"/* "$plot_png_dir"/* 2>/dev/null
# temporary working directories for gmt are within here
# gmt process create temp files which clash between process
gmt_temp="$sim_dir/gmt_wd"
mkdir -p "$gmt_temp"

# corners stored here
model_params="${vel_mod_params_dir}/model_params_nz01-h${h}"
# make temp file with simulation boundaries
if [ -e "sim.modelpath" ]; then
    \rm sim.modelpath
fi
# have to close box (copy first corner to end)
for corner in c1 c2 c3 c4 c1; do
    grep "$corner= " $model_params | gawk ' { print $2, $3 } ' >> sim.modelpath
    if [ "$corner" == "c1" ]; then
        plot_x_min=$(grep "$corner= " $model_params | gawk ' { print $2 } ')
    elif [ "$corner" == "c2" ]; then
        plot_y_max=$(grep "$corner= " $model_params | gawk ' { print $3 } ')
    elif [ "$corner" == "c3" ]; then
        plot_x_max=$(grep "$corner= " $model_params | gawk ' { print $2 } ')
    elif [ "$corner" == "c4" ]; then
        plot_y_min=$(grep "$corner= " $model_params | gawk ' { print $3 } ')
    fi
done
plot_ll_region=$plot_x_min/$plot_x_max/$plot_y_min/$plot_y_max
# create mask for simulation domain
grdmask sim.modelpath -R$plot_ll_region -I$plot_dx/$plot_dy -Gmodelmask.grd

# create color palette for plotting the topography
base_cpt=y2r_brown.cpt

# avg Lon/Lat (midpoint of the TSlice image for the geo projection)
avg_ll=(`echo $plot_x_min $plot_x_max $plot_y_min $plot_y_max | gawk '{print 0.5*($1+$2),0.5*($3+$4);}'`)
# plot projection and region in GMT format for ease of use later
att="-JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} -R${plot_ll_region}"

add_source() {
    # add the fault plane or beachball
    if [ "$plot_type" == "finitefault" ]; then
        # plot fault plane
        bash ${plot_fault_add_plane} "$1" \
                -R$plot_ll_region -JT${avg_ll[0]}/${avg_ll[1]}/${plot_x_inch} \
                $fault_file $plot_fault_line $plot_fault_top_edge $plot_fault_hyp_open
    else
        # plot beach ball (variable must be surrounded by quotes in sourced file)
        gmt psmeca -P -J -R -Sc0.15 -Ggreen -O -K << EOF >> "$1"
$plot_beachball
EOF
    fi
}

cpt=hot
min=0
extra='-I'
scale='ground velocity (cm/s)'
makecpt -C$cpt $extra -T$min/$plot_topo_a_max/$plot_topo_a_inc -A50 > $base_cpt

# set all plotting defaults to use
gmtset FONT_ANNOT_PRIMARY 16 MAP_TICK_LENGTH_PRIMARY 0.05i FONT_LABEL 16 PS_PAGE_ORIENTATION PORTRAIT MAP_FRAME_PEN 1p FORMAT_GEO_MAP D MAP_FRAME_TYPE plain FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT i

###################### BEGIN TEMPLATE ##########################
echo Creating PS Template...
plot_file_template=$gmt_temp/plot_template.ps
# specify plot and panel size (defaults 8.5 x 11)
psxy -JX0.1/0.1 -R0/0.1/0/0.1 -G255/255/255 -X0 -Y0 -K << END > "$plot_file_template"
END
# specify the X and Y offsets` for plotting (I dont really understand this yet)
psxy -V $att -L  -K -O -X$plot_x_org -Y$plot_y_org << END >> "$plot_file_template" 2>/dev/null #-W5/255/255/0
END
# simulation domain
gmt psxy -t50 sim.modelpath $att -W1.5p,black,- -L -O -K >> "$plot_file_template"
rm sim.modelpath
# include maxgrid on map
# create ground motion intensity surface from MAXGRID
if [ "$swap_bytes" -eq 1 ]; then
    xyz2grd statgrid_max.bin -Sstatgrid_max.native.bin -V -Zf 2>/dev/null
else
    cp statgrid_max.bin statgrid_max.native.bin
fi
surface statgrid_max.native.bin -Gmax.grd -I$plot_dx/$plot_dy \
        -R$plot_ll_region -T0.0 -bi3f 2>/dev/null
rm statgrid_max.native.bin
# crop to simulation domain (multiply by mask of 0 or 1, outside = 0)
grdmath max.grd modelmask.grd MUL = max.grd 2>/dev/null
# clip minimum (values below cutoff = NaN, not displayed, clear)
grdclip max.grd -Gmax.grd \
        -Sb${plot_topo_a_min}/NaN 2>/dev/null
# add resulting overlay image to plot
grdimage max.grd $att -C$base_cpt -Q -t70 -K -O >> "$plot_file_template" 2>/dev/null
# finite fault or beachball
add_source "$plot_file_template"
# have template ready and raster version
mkdir -p WEB
cp "$plot_file_template" WEB/gmt-template.ps
cp gmt.conf WEB/
cp gmt.history WEB/
psxy -V $att -L -W5,255/255/0 -O << END >>  "$plot_file_template" 2>/dev/null
END
ps2raster "$plot_file_template" -A -TG -E$plot_res -DWEB
echo Template Complete.

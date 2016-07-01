#!/usr/bin/env bash
source e3d.par

echo
echo Movie creation starts
echo
convert -delay 5 $plot_png_dir/*.png $name.gif
echo Movie creation ends


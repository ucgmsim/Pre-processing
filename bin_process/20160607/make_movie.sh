#!/usr/bin/env bash
source e3d.par
ffbuild=/nesi/home/vap30/bin/$(uname -s)/ffmpeg

# https://trac.ffmpeg.org/wiki/Create%20a%20video%20slideshow%20from%20images

# movie speed, can be fractional
# doesn't affect filesize when using qtrle codec
# 1: movie time = simulation time
# 2: movie twice as fast
speed=2
# make sure fps >= 1 (result is int)
fps=$(echo "1/($dt * $dtts) * $speed" | bc)

# create quicktime animation (flexible resolution, high quality, low size)
# alternative may be libx264 (h264 compile time external lib)
# for compatibility with old/commercial software, don't use libx265
# to scale to maximum width or hight (for variable input size):
#   -vf scale="'if(gt(a,320/240),320,-1)':'if(gt(a,320/240),-1,240)'"
#   https://trac.ffmpeg.org/wiki/Scaling%20(resizing)%20with%20ffmpeg
#   does not work with old ffmpeg versions
# -framerate: input framerate
# -i: input filename format, alternative is using glob
# -c:v: video codec, list options: 'ffmpeg -codecs | grep EV'
# -r: output framerate (keep same as input to prevent duplicated frames)
# -y: overwrite existing output without asking
echo
echo Movie creation starts
echo
$ffbuild -y -framerate $fps -i $plot_png_dir/ts-str%04d.png \
        -c:v qtrle -r $fps $name.mov 2>/dev/null
echo Movie creation ends


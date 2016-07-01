#!/bin/bash
export PATH=$PATH:/nesi/projects/nesi00213/linux_pkgs/bin
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #find the path of this file

bash $DIR/plot_ts.sh 2
bash $DIR/make_movie.sh

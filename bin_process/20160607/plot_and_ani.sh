#!/bin/bash
export PATH=$PATH:/nesi/projects/nesi00213/linux_pkgs/bin
bash $BINPROCESS/plot_ts.sh 2
bash $BINPROCESS/make_movie.sh

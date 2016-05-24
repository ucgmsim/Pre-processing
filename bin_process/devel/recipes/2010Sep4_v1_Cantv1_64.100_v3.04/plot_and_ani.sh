#!/bin/bash
export BINPROCESS=/hpc/scratch/nesi00213/Pre-processing/bin_process/20160511
export PATH=$PATH:/nesi/projects/nesi00213/linux_pkgs/bin
bash $BINPROCESS/plot_ts.sh
bash $BINPROCESS/make_movie.sh

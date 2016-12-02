#!/usr/bin/bash
source version
echo $VERSION
python /nesi/projects/nesi00213/Pre-processing/bin_process/$VERSION/submit_post_emod3d.py

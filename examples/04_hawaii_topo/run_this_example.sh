#!/bin/bash

NPROC=8

pos_str=19.5/-155.5
angle=-30

# download topography
gmt grdcut @earth_relief_01m -R-157/-152/18/21 -Ghawaii.nc

# rotate source receiver file
../../bin/surfatt_rotate_src_rec -i src_rec_file_raw.csv -a $angle -c $pos_str -o src_rec_file_rotated.csv

# rotate topography
../../bin/surfatt_rotate_topo -i hawaii.nc -a $angle -c $pos_str -x -0.75/0.8 -y -0.75/1 -o hawaii_rotated.nc

# inversion 
mpirun -np $NPROC ../../bin/surfatt_tomo -i input_params.yml

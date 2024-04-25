#!/bin/bash

NPROC=8

mkdir -p OUTPUT_FILES
cp src_rec_file_ph.csv OUTPUT_FILES/src_rec_file_forward_ph.csv
# cp src_rec_file_gr.csv OUTPUT_FILES/src_rec_file_forward_gr.csv

# create 2x3x2 checkers and forward simulate surface traveltimes
mpirun -np $NPROC ../../bin/surfatt_cb_fwd -i input_params.yml -n 2/3/2 -m 0.2 -p 0.08

# inversion 
mpirun -np $NPROC ../../bin/surfatt_tomo -i input_params.yml

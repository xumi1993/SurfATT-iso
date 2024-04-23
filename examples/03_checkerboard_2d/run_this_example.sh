#!/bin/bash

NPROC=8

# rm -f OUTPUT_FILES/
mkdir -p OUTPUT_FILES
cp src_rec_file_ph.csv OUTPUT_FILES/src_rec_file_forward_ph.csv

# inversion 
# mpirun -np $NPROC ../../bin/surfatt_tomo2d -i input_params.yml -n 2/3 -a 0.06/45 -p 0.08 -m 0.2 -f

mpirun --oversubscribe -np $NPROC ../../bin/surfatt_tomo2d -i input_params.yml
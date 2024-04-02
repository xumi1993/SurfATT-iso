#!/bin/bash

NPROC=8

# inversion 
mpirun --oversubscribe -np $NPROC ../../bin/surfatt_tomo -i input_params.yml

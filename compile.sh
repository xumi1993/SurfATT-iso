#!/bin/bash

# This script is used to build the code on the NSCC cluster
# The script is written for the GNU compilers, but can be easily modified to use the Intel compilers

# Load the required modules
module purge
module load gcc/11.2.0-nscc libfabric/1.11.0.4.125 openmpi/4.1.5-gcc11 cmake/3.23.1 hdf5/1.10.5

# Build the code
mkdir -p build && cd build && rm -rf *

# Build with Intel compilers
# CXX=icpc CC=icc FC=ifort MPIFC=mpiifort MPICC=mpiicc cmake ..

# Build with GNU compilers
CXX=g++ CC=gcc FC=gfortran MPIFC=mpifort MPICC=mpicc cmake ..

# Make the code
make clean && make -j4

module purge
#module load mpi/intel/20.0.4 hdf5/1.10.6-intel20 cmake/3.23.1 gcc/7.3.0
module load hdf5/1.10.5-icc21 mpi/intel/2021.1 cmake/3.23.1 gcc/9.3.0-new

# Build the code
mkdir -p build && cd build
rm -rf *
CXX=icpc CC=icc FC=ifort MPIFC=mpiifort MPICC=mpiicc cmake ..
make clean
make -j4

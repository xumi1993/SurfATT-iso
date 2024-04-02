# SurfATT

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)


This is an innovative package for **Surf**ace wave **A**djoint **T**ravel-time **T**omography driven by modern fortran with highlights:

- Inversion for **isotropic** or **azimuthal anisotropic** media (Hao et al., 2024b, in preparation)
- Calculation of surface wave travel time based on **Eikonal equation** with fast sweeping method ([Tong, 2021a](https://doi.org/10.1029/2021JB021818))
- Computation of isotropic and anisotropic sensitivity kernels through **adjoint method** ([Tong, 2021b](https://doi.org/10.1029/2021JB022365))
- **Multi-grid model parametrization** utilization in optimization ([Tong et al., 2019](https://doi.org/10.1093/gji/ggz151))
- Consideration of **surface topography** ([Hao et al., 2024a](https://doi.org/10.1029/2023JB027454))

## Installation

### Dependencies

- git
- Fortran (intel fortran or gfortran>=11)
- hdf5 (serial version)
- C++
- CMake (3.8 or higher)
- MPI (3.0 or higher)
- fypp (required by external lib [fortran-stdlib](https://github.com/fortran-lang/stdlib))

### Install dependencies

#### On local computer

- For Ununtu/Debian users

```
sudo apt install build-essential gfortran cmake libopenmpi-dev openmpi-bin libhdf5-dev
```

- For Centos/Fedora users

```
sudo dnf install gcc gcc-gfortran gcc-c++ cmake openmpi-devel hdf5-devel
```

- For MacOS users

```
brew install gcc cmake open-mpi hdf5
```

#### On HPC

The `module-environment` is extensively utilized for managing modules that are requisite for SurfATT. Dependencies can be loaded as illustrated below

```
module load mpi/intel/20.0.4 hdf5/1.10.6-intel20 cmake/2.23.1 
```

- In the case of HPCs where `hdf5` for Fortran might not be pre-installed, it is possible to manually compile and install `hdf5`.

    ``` bash
    cd ./SurfATT/external_libs

    # make a local install pass
    mkdir local_hdf5 

    # download hdf5 source
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.3/src/hdf5-1.13.3.tar.gz
    
    # extract the downloaded directory
    tar -xvf hdf5-1.13.3.tar.gz
    cd hdf5-1.13.3/
    
    # configure
    CC=mpicc CXX=mpic++ FC=gfortran ./configure --enable-unsupported --enable-shared --enable-fortran --prefix=$(pwd)/../local_hdf5
    
    # make
    make -j && make install
    ```

#### Install `fypp`

`fypp` is required by [fortran-stdlib](https://github.com/fortran-lang/stdlib), which can be easily installed via the `pip`:

```
pip install fypp
```

## Download and Compile SurfATT

### Clone the source codes

```
git clone --recursive https://github.com/xumi1993/SurfATT
```
`--recursive` is required in this command, because dependence of `yaml-cpp`, `h5fortran`, and `fortran-stdlib` will be download automatically.

### Compile source codes with cmake

```
mkdir build && cd build
cmake .. && make -j4
```

For MacOS users or if you have manually installed compilers, you can specify path to compile by add following options before cmake:

```
CC=gcc-12 CXX=g++-12 FC=gfortran-12 MPIFC=mpif90 cmake ..
make -j4
```

For users who opt for a manual installation of `hdf5`, the local `hdf5` directory should be specified before running `cmake`.
```
CC=gcc CXX=g++ FC=gfortran MPIFC=mpif90 HDF5_ROOT=$(pwd)/../external_libs/local_hdf5  cmake ..
make -j4
```

## How to use SurfATT
The executable file `bin/surfatt_tomo` for inverting surface dispersion data for S-wave velocity can be run with `mpirun` as:
```
mpirun -np 4 bin/attsurf_tomo -i input_params.yml
```

### A quick example
A case named `test/00_checkerboard_iso` presents an example of inversion for 2x3x2 checkers using ambient noise surface wave data from 25 stations. execute run_this_example.sh to run this example under 5 processors.
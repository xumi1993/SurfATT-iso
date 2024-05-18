# SurfATT

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![License](https://img.shields.io/github/license/xumi1993/seispy)]()

This is an innovative package for **Surf**ace wave **A**djoint **T**ravel-time **T**omography driven by modern fortran with highlights:

- Calculation of surface wave travel time based on **Eikonal equation** with fast sweeping method ([Tong, 2021a](https://doi.org/10.1029/2021JB021818))
- Computation of sensitivity kernels through **adjoint method** ([Tong, 2021b](https://doi.org/10.1029/2021JB022365))
- **Multi-grid model parametrization** utilization in optimization ([Tong et al., 2019](https://doi.org/10.1093/gji/ggz151))
- Consideration of **surface topography** ([Hao et al., 2024a](https://doi.org/10.1029/2023JB027454))
![Fig2](https://github.com/xumi1993/SurfATT-iso/assets/7437523/f9a0155b-7b83-4970-914d-f13dc42b11e5)

## Installation

Please refer to the [installation guide](https://surfatt.xumijian.me/installation/dependence.html) for detailed instructions.

## How to use SurfATT
The executable file `bin/surfatt_tomo` for inverting surface dispersion data for S-wave velocity can be run with `mpirun` as:
```
mpirun -np 4 bin/surfatt_tomo -i input_params.yml
```

### A quick example
A case named `test/00_checkerboard_iso` presents an example of inversion for 2x3x2 checkers using ambient noise surface wave data from 25 stations. execute `run_this_example.sh` to run this example under 8 processors.

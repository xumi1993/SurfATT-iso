# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

add_library(fortran-yaml OBJECT yaml_types.f90 yaml.f90 yaml.cpp)
target_link_libraries(fortran-yaml yaml-cpp)
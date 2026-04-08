# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

add_library(fortran-yaml OBJECT yaml_types.f90 yaml.f90 yaml.cpp)
# Use YAML_CPP_TARGET set by parent CMakeLists.txt, or determine if not set
if (NOT DEFINED YAML_CPP_TARGET)
  # This shouldn't happen if parent is properly configured
  message(FATAL_ERROR "YAML_CPP_TARGET not defined by parent CMakeLists.txt")
endif()
target_link_libraries(fortran-yaml ${YAML_CPP_TARGET})
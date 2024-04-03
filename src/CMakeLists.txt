# Shared library
add_subdirectory("${SOURCE_DIR}/shared")

# fortran-yaml-cpp
add_subdirectory("${SOURCE_DIR}/fortran-yaml")

# fortran-csv
add_subdirectory("${SOURCE_DIR}/fortran-csv")

# surfker
add_subdirectory("${SOURCE_DIR}/surfker")

file(GLOB SOURCES "${SOURCE_DIR}/*.f90")
file(GLOB interp_sources "${SOURCE_DIR}/interp2d_scatter/*.f90")
list(APPEND SOURCES ${interp_sources})
add_library(surfatt_modules OBJECT ${SOURCES})
target_link_libraries(surfatt_modules PRIVATE h5fortran fortran-yaml fortran-csv shared_obj surfker_obj)

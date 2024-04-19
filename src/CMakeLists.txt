# Shared library
add_subdirectory("${SOURCE_DIR}/shared")

# fortran-yaml-cpp
add_subdirectory("${SOURCE_DIR}/fortran-yaml")

# fortran-csv
add_subdirectory("${SOURCE_DIR}/fortran-csv")

# surfker
add_subdirectory("${SOURCE_DIR}/surfker")

file(GLOB SOURCES "${SOURCE_DIR}/*.f90")

set(FILTERED_SOURCE_FILES "")
foreach(FILE ${SOURCES})
    if (NOT ${FILE} MATCHES "surfatt.*")
        list(APPEND FILTERED_SOURCE_FILES ${FILE})
    endif()
endforeach()

add_library(surfatt_modules OBJECT ${FILTERED_SOURCE_FILES})
# target_link_libraries(surfatt_modules PRIVATE fortran-yaml fortran-csv shared_obj surfker_obj ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
target_link_libraries(surfatt_modules PRIVATE fortran-yaml fortran-csv shared_obj surfker_obj ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} MPI::MPI_Fortran)

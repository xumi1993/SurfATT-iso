cmake_minimum_required(VERSION 3.8)

# set the project name
set(CMAKE_PROJECT_NAME "SurfATT")
project(${CMAKE_PROJECT_NAME} VERSION 1.4.0 LANGUAGES Fortran CXX)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(MPI_Fortran_COMPILE_FLAGS "-ffixed-line-length-none -O3")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(MPI_Fortran_COMPILE_FLAGS "-132 -xHost -O3")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")

set(CMAKE_CXX_STANDARD 17)
set(BUILD_TESTING OFF)
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# setup mpi
find_package(MPI REQUIRED)
message(STATUS "MPI_FOUND: ${MPI_FOUND}")
message(STATUS "MPI_VERSION: ${MPI_VERSION}")

find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
include_directories(${HDF5_INCLUDE_DIRS} ${HDF5_HL_INCLUDE_DIR})
message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
message(STATUS "HDF5_VERSION: ${HDF5_VERSION}")

# add sub modules
add_subdirectory("${PROJECT_SOURCE_DIR}/external_libs")

# surfatt_modules
add_subdirectory("${SOURCE_DIR}")

# setup include directory
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

set(APP_SOURCES
${SOURCE_DIR}/surfatt_tomo.f90
${SOURCE_DIR}/surfatt_cb_fwd.f90
${SOURCE_DIR}/surfatt_tomo2d.f90
${SOURCE_DIR}/surfatt_rotate_src_rec.f90
${SOURCE_DIR}/surfatt_rotate_topo.f90
  )
foreach(execsourcefile ${APP_SOURCES})
  # get exec name
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  # add executable
  add_executable(${EXEC_NAME} ${execsourcefile})

  # link MPI
  target_link_libraries(${EXEC_NAME} PUBLIC MPI::MPI_Fortran)
  
  # link yaml-app:
  target_link_libraries(${EXEC_NAME} PUBLIC fortran-yaml)

  # link fortran-csv
  target_link_libraries(${EXEC_NAME} PUBLIC fortran-csv)
  
  # link h5fortran
  target_link_libraries(${EXEC_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

  # link surfatt_modules
  target_link_libraries(${EXEC_NAME} PUBLIC shared_obj surfker_obj surfatt_modules)

  # set linker language
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

endforeach()

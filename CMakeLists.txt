cmake_minimum_required(VERSION 3.8)

# set the project name
set(CMAKE_PROJECT_NAME "SurfATT")
project(${CMAKE_PROJECT_NAME} VERSION 1.4.0 LANGUAGES CXX Fortran)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(MPI_Fortran_COMPILE_FLAGS "-ffixed-line-length-none -fmax-errors=10 -ftree-vectorize -march=native -O3 -ffast-math")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  set(MPI_Fortran_COMPILE_FLAGS "-ffixed-line-length-none -xHost -O3")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")

set(CMAKE_CXX_STANDARD 17)
set(BUILD_TESTING OFF)
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

# setup mpi
find_package(MPI REQUIRED)
message(STATUS "MPI_FOUND: ${MPI_FOUND}")
message(STATUS "MPI_VERSION: ${MPI_VERSION}")

# add sub modules
add_subdirectory("${PROJECT_SOURCE_DIR}/external_libs")

# surfatt_modules
add_subdirectory("${SOURCE_DIR}")

# setup include directory
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

set(APP_SOURCES
  src/surfatt_tomo.f95
  src/surfatt_cb_fwd.f95
  src/surfatt_tomo2d.f95
  src/surfatt_rotate_src_rec.f95
  src/surfatt_rotate_topo.f95
  )
foreach(execsourcefile ${APP_SOURCES})
  # get exec name
  get_filename_component(EXEC_NAME ${execsourcefile} NAME_WE)

  # add executable
  add_executable(${EXEC_NAME} ${execsourcefile})

  # link MPI
  target_link_libraries(${EXEC_NAME} PUBLIC MPI::MPI_Fortran)
  
  # link fortran stdlib
  # target_link_libraries(${EXEC_NAME} PUBLIC fortran_stdlib)

  # link yaml-app:
  target_link_libraries(${EXEC_NAME} PUBLIC fortran-yaml)

  # link fortran-csv
  target_link_libraries(${EXEC_NAME} PUBLIC fortran-csv)
  
  # link h5fortran
  target_link_libraries(${EXEC_NAME} PUBLIC h5fortran)

  # link surfatt_modules
  target_link_libraries(${EXEC_NAME} PUBLIC shared_obj surfker_obj surfatt_modules)

endforeach()

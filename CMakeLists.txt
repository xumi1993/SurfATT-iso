cmake_minimum_required(VERSION 3.10)

# set the project name
set(CMAKE_PROJECT_NAME "SurfATT")
project(${CMAKE_PROJECT_NAME} VERSION 1.6.3 LANGUAGES Fortran CXX)

# set the version number and git commit hash
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
  ${PROJECT_SOURCE_DIR}/include/version.h.in
  ${PROJECT_SOURCE_DIR}/include/version.h
)

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
find_package(yaml-cpp 0.8 QUIET)
if (yaml-cpp_FOUND)
  message(STATUS "yaml-cpp found via find_package")
  # Determine the correct target name (newer versions use yaml-cpp::yaml-cpp)
  if (TARGET yaml-cpp::yaml-cpp)
    set(YAML_CPP_TARGET yaml-cpp::yaml-cpp)
    message(STATUS "Using target: yaml-cpp::yaml-cpp")
  elseif (TARGET yaml-cpp)
    set(YAML_CPP_TARGET yaml-cpp)
    message(STATUS "Using target: yaml-cpp")
  else()
    message(FATAL_ERROR "yaml-cpp found but no target available")
  endif()
else()
  message(STATUS "yaml-cpp not found via find_package. Trying pkg-config...")
  find_package(PkgConfig QUIET)
  if (PkgConfig_FOUND)
    pkg_check_modules(YAML_CPP yaml-cpp QUIET)
    if (YAML_CPP_FOUND)
      message(STATUS "yaml-cpp found via pkg-config")
      message(STATUS "YAML_CPP_LIBRARIES: ${YAML_CPP_LIBRARIES}")
      message(STATUS "YAML_CPP_LIBRARY_DIRS: ${YAML_CPP_LIBRARY_DIRS}")
      # Find the actual library file
      find_library(YAML_CPP_LIB NAMES yaml-cpp PATHS ${YAML_CPP_LIBRARY_DIRS})
      if (YAML_CPP_LIB)
        message(STATUS "Found yaml-cpp library: ${YAML_CPP_LIB}")
        # Create an imported library target from pkg-config info
        add_library(yaml-cpp-pc SHARED IMPORTED)
        set_target_properties(yaml-cpp-pc PROPERTIES
          IMPORTED_LOCATION "${YAML_CPP_LIB}"
          INTERFACE_INCLUDE_DIRECTORIES "${YAML_CPP_INCLUDE_DIRS}"
        )
        set(YAML_CPP_TARGET yaml-cpp-pc)
      else()
        message(FATAL_ERROR "yaml-cpp pkg-config found but library not located")
      endif()
    else()
      message(STATUS "yaml-cpp not found via pkg-config either. Using external_libs/yaml-cpp...")
      add_subdirectory("${PROJECT_SOURCE_DIR}/external_libs")
      set(YAML_CPP_TARGET yaml-cpp)
    endif()
  else()
    message(STATUS "PkgConfig not available. Using external_libs/yaml-cpp...")
    add_subdirectory("${PROJECT_SOURCE_DIR}/external_libs")
    set(YAML_CPP_TARGET yaml-cpp)
  endif()
endif()

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
${SOURCE_DIR}/surfatt_rotate_model.f90
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
  
  # link yaml-cpp library explicitly (to resolve symbols from yaml.cpp.o)
  target_link_libraries(${EXEC_NAME} PUBLIC ${YAML_CPP_TARGET})

  # link fortran-csv
  target_link_libraries(${EXEC_NAME} PUBLIC fortran-csv)
  
  # link h5fortran
  target_link_libraries(${EXEC_NAME} PUBLIC ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

  # link surfatt_modules
  target_link_libraries(${EXEC_NAME} PUBLIC shared_obj surfker_obj surfatt_modules)

  # set linker language
  set_target_properties(${EXEC_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

endforeach()

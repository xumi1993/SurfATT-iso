# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

file(GLOB SHARED_SOURCES "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(shared_obj OBJECT ${SHARED_SOURCES})
target_link_libraries(shared_obj PRIVATE MPI::MPI_Fortran fortran-csv)
# target_link_libraries(fortran-yaml)
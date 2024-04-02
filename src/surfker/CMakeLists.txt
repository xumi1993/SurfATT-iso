# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

file(GLOB SURFKER_SOURCES "*.f*")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(surfker_obj OBJECT ${SURFKER_SOURCES})
target_link_libraries(surfker_obj PRIVATE shared_obj)

# target_link_libraries(fortran-yaml)
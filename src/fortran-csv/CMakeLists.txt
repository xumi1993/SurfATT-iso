# include(FortranCInterface)
# FortranCInterface_VERIFY(CXX)

set(FCSV_SOURCES
    "csv_kinds.f90"
    "csv_parameters.f90"
    "csv_utilities.f90"
    "csv_module.F90"
)

add_library(fortran-csv OBJECT ${FCSV_SOURCES})
# target_link_libraries(fortran-yaml)
# Once done, this will define
#  P4EST_FOUND - true if p4est has been found
#  P4EST_INCLUDE_DIR - the p4est include dir
#  P4EST_LIBRARIES - names of p4est libraries
#  P4EST_LINK_DIRECTORY - location of p4est libraries

find_path(P4EST_INCLUDE_DIR p4est.h HINTS ${P4EST_ROOT}/../include PATH_SUFFIXES include)

find_library(P4EST_LIBRARIES NAMES p4est HINTS ${P4EST_DIR}/../lib PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(P4EST DEFAULT_MSG P4EST_INCLUDE_DIR P4EST_LIBRARIES)

if(NOT SC_LIBRARIES)
    message(FATAL_ERROR "p4est libraries not found")
endif()

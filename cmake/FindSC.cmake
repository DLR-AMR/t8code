# Once done, this will define
#  SC_FOUND - true if SC has been found
#  SC_INCLUDE_DIR - the SC include dir
#  SC_LIBRARIES - names of SC libraries
#  SC_LINK_DIRECTORY - location of SC libraries

find_path(SC_INCLUDE_DIR sc.h HINTS ${SC_DIR}/../include PATH_SUFFIXES include)

find_library(SC_LIBRARIES NAMES sc HINTS ${SC_DIR}/../lib PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SC DEFAULT_MSG SC_INCLUDE_DIR SC_LIBRARIES)

if(NOT SC_LIBRARIES)
    message(FATAL_ERROR "SC libraries not found")
endif()

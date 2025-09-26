# This file is inspired by the project cmake-modules (https://github.com/bilke/cmake-modules) 
# of Lars Bilke but heavily adapted and rewritten!
# The project is licensed with the BSD license. 
# Copyright (c) 2012 - 2017, Lars Bilke
# All rights reserved.
# See also the license information in CodeCoverage.license.

# This file defines and calls functions to collect code coverage information.
# A code coverage report is created in a folder "coverage".
# Lcov and genhtml are used to generate the report and should be available.

#
# Adds compiler flags necessary to be able to collect coverage information.
# 
function(append_coverage_compiler_flags)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g --coverage -O0" PARENT_SCOPE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g --coverage -O0" PARENT_SCOPE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage" PARENT_SCOPE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage" PARENT_SCOPE)
endfunction() # append_coverage_compiler_flags

#
# Defines a target to collect code coverage information.
# Lcov and genhtml are used to generate a coverage report and should be available.
#
# \param NAME: New target name.
# \param EXCLUDE: "src/dir1/*" "src/dir2/*": Folders or files to exclude from coverage report.
# \param LCOV_ARGS: lcov arguments. 
# 
function(setup_target_for_coverage)
  set(options "")
  set(oneValueArgs NAME)
  set(multiValueArgs EXCLUDE LCOV_ARGS)
  cmake_parse_arguments(Coverage "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Find lcov and genhtml and mark as advanced options.
  find_program( LCOV_PATH  NAMES lcov lcov.bat lcov.exe lcov.perl)
  find_program( GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat )
  mark_as_advanced(
    LCOV_PATH 
    GENHTML_PATH
  )
  if(NOT LCOV_PATH)
    message(FATAL_ERROR "LCOV not found but required for code coverage report! Aborting...")
  endif()
  if(NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found but required for code coverage report! Aborting...")
  endif()

  # Project dir is used as base dir for lcov.
  set(PROJECT_DIR ${PROJECT_SOURCE_DIR})

  # If possible, use more than one job to execute the test suite.
  cmake_host_system_information(
    RESULT N
    QUERY NUMBER_OF_LOGICAL_CORES
  )
  # If the call leads to an error, use just one job.
  if(N EQUAL 0)
    set(N 1)
  endif()

  # Setup coverage target.
  add_custom_target(${Coverage_NAME}
    # Cleanup lcov from prior usages.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${PROJECT_DIR} --zerocounters
    # Create baseline to make sure untouched files show up in the report.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${PROJECT_DIR} --capture --initial -output-file ${Coverage_NAME}.base

    # Run tests and collect coverage information.
    COMMAND ctest -T Test -T Coverage -j ${N}
    
    # Generate report using lcov.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${PROJECT_DIR} --capture --output-file ${Coverage_NAME}.capture
    # Add baseline counters created above.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --add-tracefile ${Coverage_NAME}.base --add-tracefile ${Coverage_NAME}.capture --output-file ${Coverage_NAME}.total
    # Filter collected data using the variable Coverage_EXCLUDE and generate final coverage report.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --remove ${Coverage_NAME}.total ${Coverage_EXCLUDE} --output-file ${Coverage_NAME}.info

    # Generate HTML file using genhtml.
    COMMAND ${GENHTML_PATH} -o ${Coverage_NAME} ${Coverage_NAME}.info

    BYPRODUCTS
      ${Coverage_NAME}.base
      ${Coverage_NAME}.capture
      ${Coverage_NAME}.total
      ${Coverage_NAME}.info
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    VERBATIM
    COMMENT "Resetting code coverage counters and generate new coverage report."
    )

endfunction() # setup_target_for_coverage

# Call functions defined above with customized arguments.
append_coverage_compiler_flags()
setup_target_for_coverage(
  NAME coverage
  EXCLUDE "${CMAKE_SOURCE_DIR}/sc*" "${CMAKE_SOURCE_DIR}/p4est*" "${CMAKE_SOURCE_DIR}/test*" "${CMAKE_SOURCE_DIR}/thirdparty*" "${CMAKE_SOURCE_DIR}/tutorials*" "${CMAKE_SOURCE_DIR}/example*" "${CMAKE_SOURCE_DIR}/benchmarks*"
  LCOV_ARGS --no-external --ignore-errors gcov,mismatch --demangle-cpp
)
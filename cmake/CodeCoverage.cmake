function(append_coverage_compiler_flags)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g --coverage" PARENT_SCOPE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g --coverage" PARENT_SCOPE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage" PARENT_SCOPE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage" PARENT_SCOPE)
endfunction() # append_coverage_compiler_flags

#
# Defines a target to collect code coverage information.
# Lcov and genhtml are used to generate a coverage report and should be available.
#
# \param NAME: New target name.
# \param EXCLUDE: "src/dir1/*" "src/dir2/*": Folder or Files to exclude from coverage report.
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
    message(FATAL_ERROR "LCOV not found! Aborting...")
  endif()
  if(NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found! Aborting...")
  endif()

  # Project dir is used as based dir for lcov.
  set(PROJECT_DIR ${PROJECT_SOURCE_DIR})

  # If possible, use more than one job to execute the test suite.
  cmake_host_system_information(
    RESULT N
    QUERY NUMBER_OF_PHYSICAL_CORES
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
    COMMAND ctest -T Test -T Coverage -j 4
    
    # Generate report using lcov
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${PROJECT_DIR} --capture --output-file ${Coverage_NAME}.capture
    # Add baseline counters created above.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --add-tracefile ${Coverage_NAME}.base --add-tracefile ${Coverage_NAME}.capture --output-file ${Coverage_NAME}.total
    # Filter collected data using the variable Coverage_EXCLUDE and generate final coverage report.
    COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --remove ${Coverage_NAME}.total ${Coverage_EXCLUDE} --output-file ${Coverage_NAME}.info

    # Generate HTML file using genhtml.
    COMMAND ${GENHTML_PATH} -o ${Coverage_NAME} ${Coverage_NAME}.info
    # Show result in terminal.
    COMMAND ${LCOV_PATH} --list ${Coverage_NAME}.info

    # Set output files as GENERATED
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
  EXCLUDE "${CMAKE_SOURCE_DIR}/sc*"
  LCOV_ARGS "--no-external"
)
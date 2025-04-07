function(append_coverage_compiler_flags)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g --coverage" PARENT_SCOPE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g --coverage" PARENT_SCOPE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage" PARENT_SCOPE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage" PARENT_SCOPE)
endfunction() # append_coverage_compiler_flags

# Defines a target for running and collection code coverage information.
#
# define_coverage_target(
#     NAME testrunner_coverage                    # New target name
#     BASE_DIRECTORY "../"                        # Base directory for report
#                                                 #  (defaults to PROJECT_SOURCE_DIR)
#     EXCLUDE "src/dir1/*" "src/dir2/*"           # Patterns to exclude (can be relative
#                                                 #  to BASE_DIRECTORY, with CMake 3.4+)
# )
function(setup_target_for_coverage)
  find_program( LCOV_PATH  NAMES lcov lcov.bat lcov.exe lcov.perl)
  find_program( GENHTML_PATH NAMES genhtml genhtml.perl genhtml.bat )

  mark_as_advanced(
    LCOV_PATH 
    GENHTML_PATH)

  if(NOT LCOV_PATH)
    message(FATAL_ERROR "LCOV not found! Aborting...")
  endif()

  if(NOT GENHTML_PATH)
    message(FATAL_ERROR "genhtml not found! Aborting...")
  endif()

    # Collect excludes (CMake 3.4+: Also compute absolute paths)
    # set(LCOV_EXCLUDES "")
    # foreach(EXCLUDE ${Coverage_EXCLUDE} ${COVERAGE_EXCLUDES} ${COVERAGE_LCOV_EXCLUDES})
    #     if(CMAKE_VERSION VERSION_GREATER 3.4)
    #         get_filename_component(EXCLUDE ${EXCLUDE} ABSOLUTE BASE_DIR ${BASEDIR})
    #     endif()
    #     list(APPEND LCOV_EXCLUDES "${EXCLUDE}")
    # endforeach()
    # list(REMOVE_DUPLICATES LCOV_EXCLUDES)

    # message(${LCOV_EXCLUDES})
    set(BASEDIR ${PROJECT_SOURCE_DIR})
    set(Coverage_NAME coverage)
    set(LCOV_EXCLUDES "${CMAKE_SOURCE_DIR}/sc*")
    set(Coverage_LCOV_ARGS "--no-external")

    cmake_host_system_information(RESULT N
                              QUERY NUMBER_OF_PHYSICAL_CORES)
                              if(N EQUAL 0)
                              set(N 1)
                              endif()

    # Setup target
    add_custom_target(${Coverage_NAME}

        # Cleanup lcov
        COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${BASEDIR} --zerocounters
        # Create baseline to make sure untouched files show up in the report
        COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${BASEDIR} --capture --initial -output-file ${Coverage_NAME}.base

        # Run tests
    
        COMMAND ctest -T Test -T Coverage -j ${N}
        

        # Capturing lcov counters and generating report
        COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --directory . --base-directory ${BASEDIR} --capture --output-file ${Coverage_NAME}.capture
        # add baseline counters
        COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --add-tracefile ${Coverage_NAME}.base --add-tracefile ${Coverage_NAME}.capture --output-file ${Coverage_NAME}.total
        # filter collected data to final coverage report
        COMMAND ${LCOV_PATH} ${Coverage_LCOV_ARGS} --remove ${Coverage_NAME}.total ${LCOV_EXCLUDES} --output-file ${Coverage_NAME}.info

        # Generate HTML output
        COMMAND ${GENHTML_PATH} -o ${Coverage_NAME} ${Coverage_NAME}.info

        # Show in terminal 
        COMMAND ${LCOV_PATH} --list ${Coverage_NAME}.info

        # Set output files as GENERATED
        BYPRODUCTS
            ${Coverage_NAME}.base
            ${Coverage_NAME}.capture
            ${Coverage_NAME}.total
            ${Coverage_NAME}.info
            # ${Coverage_NAME}  # report directory

        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        VERBATIM
        COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
    )

endfunction() # setup_target_for_coverage

append_coverage_compiler_flags()
setup_target_for_coverage()
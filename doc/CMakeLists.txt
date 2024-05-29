find_package(Doxygen REQUIRED)

if (DOXYGEN_FOUND)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
    set(DOXYFILE_IN ${PROJECT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYFILE_OUT ${PROJECT_BINARY_DIR}/doc/Doxyfile)

    #The following options are used inside the do
    set(top_srcdir ${PROJECT_SOURCE_DIR})
    set(top_builddir ${PROJECT_BINARY_DIR})

    set(PACKAGE_NAME ${PROJECT_NAME})
    set(VERSION ${T8_VERSION})


    #Replace variables inside @@ with the current values
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    add_custom_command( OUTPUT ${DOXYGEN_INDEX_FILE}
                        DEPENDS ${T8_PUBLIC_HEADERS}
                        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT} 
                        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN} 
                        COMMENT "Generating documentation")

    add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
else (DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
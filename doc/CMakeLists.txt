find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

if (SPHINX_FOUND)
    set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/doc/source)
    set(SPHINX_BUILD ${PROJECT_BINARY_DIR}/doc/sphinx)

    add_custom_target(Sphinx ALL
    COMMAND ${SPHINX_EXECUTABLE} -b html
    # Tell Breathe where to find the Doxygen output
    -Dbreathe_projects.CatCutifier=${DOXYGEN_OUTPUT_DIR}
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")

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
else (SPHINX_FOUND)
    message(FATAL_ERROR "Doxygen need to be installed to generate the doxygen documentation")
endif (SPHINX_FOUND)
#  This file is part of t8code.
#  t8code is a C library to manage a collection (a forest) of multiple
#  connected adaptive space-trees of general element types in parallel.
#
#  Copyright (C) 2025 the developers
#
#  t8code is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  t8code is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with t8code; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

if(T8CODE_BUILD_DOCUMENTATION)
    find_package(Doxygen REQUIRED)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc)
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
                            COMMENT "Generating documentation"
                            VERBATIM)

        add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

    else (DOXYGEN_FOUND)
        message(FATAL_ERROR "Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
endif(T8CODE_BUILD_DOCUMENTATION)

if (T8CODE_BUILD_DOCUMENTATION_SPHINX)
    find_package(Sphinx REQUIRED)

    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    if (SPHINX_FOUND)
        set(SPHINX_SOURCE ${PROJECT_SOURCE_DIR}/doc/source)
        set(SPHINX_BUILD ${PROJECT_BINARY_DIR}/doc/sphinx)
        set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

        set(CONF_FILE_IN ${SPHINX_SOURCE}/conf.in)
        set(CONF_FILE_OUT ${SPHINX_SOURCE}/conf.py)

        set(top_srcdir ${PROJECT_SOURCE_DIR})

        configure_file(${CONF_FILE_IN} ${CONF_FILE_OUT} @ONLY)

        # Only regenerate Sphinx when:
        # - Doxygen has rerun
        # - Our doc files have been updated
        # - The Sphinx config has been updated
        add_custom_command( OUTPUT ${SPHINX_INDEX_FILE}
                            COMMAND ${SPHINX_EXECUTABLE} -b html
                            # Tell Breathe where to find the Doxygen output
                            -Dbreathe_projects.T8code=${DOXYGEN_OUTPUT_DIR}/xml
                            ${SPHINX_SOURCE} ${SPHINX_BUILD}
                            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                            DEPENDS
                            # Other docs files you want to track should go here (or in some variable)
                            ${SPHINX_SOURCE}/index.rst
                            ${DOXYGEN_INDEX_FILE}
                            MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                            COMMENT "Generating documentation with Sphinx")


        # Nice named target so we can run the job easily
        add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})
    else (SPHINX_FOUND)
        message(FATAL_ERROR "Sphinx need to be installed to generate the doxygen documentation")
    endif (SPHINX_FOUND)
endif(T8CODE_BUILD_DOCUMENTATION_SPHINX)

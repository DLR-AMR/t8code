find_package( Git REQUIRED )

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    # See `scr/t8_version.h` for the documentation of following definitions.

    execute_process( COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
                    COMMAND cut -c 2-
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                    OUTPUT_VARIABLE T8CODE_VERSION_RAW
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    execute_process( COMMAND echo ${T8CODE_VERSION_RAW}
                    COMMAND cut -d- -f1
                    OUTPUT_VARIABLE T8CODE_VERSION_NUMBERS
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    execute_process( COMMAND echo ${T8CODE_VERSION_RAW}
                COMMAND cut -d- -f2-
                OUTPUT_VARIABLE T8CODE_VERSION_POINT
                OUTPUT_STRIP_TRAILING_WHITESPACE )

    # To reuse the version in other CMakeLists.
    else()
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" VERSION_CONTENT)
        # Extract the version number 
        string(REGEX MATCH "Version ([0-9])\.([0-9]+)\.([0-9]+)" VERSION_MATCH "${VERSION_CONTENT}" )
        if (VERSION_MATCH)
            # The version number will be in ${CMAKE_MATCH_1}
            set(T8CODE_VERSION_NUMBERS "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
            set(T8CODE_VERSION_RAW "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
            message(STATUS "Extracted Version: ${T8CODE_VERSION_NUMBERS}")
        else()
            message(WARNING "Version number not found in version.txt")
        endif()
    else()
        message(WARNING "Version information not found")
    endif()
endif()


execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f1
                OUTPUT_VARIABLE T8CODE_VERSION_MAJOR
                OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f2
                OUTPUT_VARIABLE T8CODE_VERSION_MINOR
                OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_NUMBERS}
                COMMAND cut -d. -f3
                OUTPUT_VARIABLE T8CODE_VERSION_PATCH
                OUTPUT_STRIP_TRAILING_WHITESPACE )
set(T8_VERSION ${T8CODE_VERSION_NUMBERS} CACHE INTERNAL "")



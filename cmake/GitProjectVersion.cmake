find_package( Git REQUIRED )
execute_process( COMMAND ${GIT_EXECUTABLE} describe --tags
                 COMMAND cut -c 2-
                 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                 OUTPUT_VARIABLE T8CODE_VERSION_RAW
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_RAW}
                 COMMAND cut -d- -f1
                 OUTPUT_VARIABLE T8CODE_VERSION
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION_RAW}
                 COMMAND cut -d- -f2-
                 OUTPUT_VARIABLE T8CODE_VERSION_COMMITS
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION}
                 COMMAND cut -d. -f1
                 OUTPUT_VARIABLE T8CODE_VERSION_MAJOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION}
                 COMMAND cut -d. -f2
                 OUTPUT_VARIABLE T8CODE_VERSION_MINOR
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND echo ${T8CODE_VERSION}
                 COMMAND cut -d. -f3
                 OUTPUT_VARIABLE T8CODE_VERSION_PATCH
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

# To reuse the version in other CMakeLists.
set(T8_VERSION ${T8CODE_VERSION} CACHE INTERNAL "")

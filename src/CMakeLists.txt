file(GLOB T8_SRC *.c *.cxx *.h *.hxx)
file(GLOB T8_CMESH_SRC t8_cmesh/*.c t8_cmesh/*.cxx t8_cmesh/*.h t8_cmesh/*.hxx)
file(GLOB T8_DATA_SRC t8_data/*.c t8_data/*.cxx t8_data/*.h t8_data/*.hxx)
file(GLOB T8_FOREST_SRC t8_forest/*.c t8_forest/*.cxx t8_forest/*.h t8_forest*.hxx)
file(GLOB T8_GEOMETRY_SRC t8_geometry/*.c t8_geometry/*.cxx t8_geometry/*.h t8_geometry/*.hxx
        t8_geometry/t8_geometry_implementations/*.c t8_geometry/t8_geometry_implementations/*.cxx
        t8_geometry/t8_geometry_implementations/*.h t8_geometry/t8_geometry_implementations/*.hxx)
file(GLOB T8_SCHEMES_SRC t8_schemes/t8_default/*.c t8_schemes/t8_default/*.cxx t8_schemes/t8_default/*.h t8_schemes/t8_default/*.hxx
        t8_schemes/t8_default/t8_default_common/*.c t8_schemes/t8_default/t8_default_common/*.cxx
	t8_schemes/t8_default/t8_default_common/*.h t8_schemes/t8_default/t8_default_common/*.hxx
        t8_schemes/t8_default/t8_default_hex/*.c t8_schemes/t8_default/t8_default_hex/*.cxx
	t8_schemes/t8_default/t8_default_hex/*.h t8_schemes/t8_default/t8_default_hex/*.hxx
        t8_schemes/t8_default/t8_default_line/*.c t8_schemes/t8_default/t8_default_line/*.cxx
	t8_schemes/t8_default/t8_default_line/*.h t8_schemes/t8_default/t8_default_line/*.hxx
        t8_schemes/t8_default/t8_default_prism/*.c t8_schemes/t8_default/t8_default_prism/*.cxx
	t8_schemes/t8_default/t8_default_prism/*.h t8_schemes/t8_default/t8_default_prism/*.hxx
        t8_schemes/t8_default/t8_default_pyramid/*.c t8_schemes/t8_default/t8_default_pyramid/*.cxx
	t8_schemes/t8_default/t8_default_pyramid/*.h t8_schemes/t8_default/t8_default_pyramid/*.hxx
        t8_schemes/t8_default/t8_default_quad/*.c t8_schemes/t8_default/t8_default_quad/*.cxx
	t8_schemes/t8_default/t8_default_quad/*.h t8_schemes/t8_default/t8_default_quad/*.hxx
        t8_schemes/t8_default/t8_default_tet/*.c t8_schemes/t8_default/t8_default_tet/*.cxx
	t8_schemes/t8_default/t8_default_tet/*.h t8_schemes/t8_default/t8_default_tet/*.hxx
        t8_schemes/t8_default/t8_default_tri/*.c t8_schemes/t8_default/t8_default_tri/*.cxx
	t8_schemes/t8_default/t8_default_tri/*.h t8_schemes/t8_default/t8_default_tri/*.hxx
        t8_schemes/t8_default/t8_default_vertex/*.c t8_schemes/t8_default/t8_default_vertex/*.cxx
	t8_schemes/t8_default/t8_default_vertex/*.h t8_schemes/t8_default/t8_default_vertex/*.hxx)
file(GLOB T8_VTK_SRC t8_vtk/*.c t8_vtk/*.cxx t8_vtk/*.h t8_vtk/*.hxx)

if ( ${T8CODE_BUILD_AS_SHARED_LIBRARY} )
    add_library( T8 SHARED
        ${T8_SRC} ${T8_CMESH_SRC} ${T8_DATA_SRC} ${T8_FOREST_SRC} ${T8_GEOMETRY_SRC} ${T8_SCHEMES_SRC} ${T8_VTK_SRC}
    )
    set_target_properties( T8 PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
    add_library( T8 STATIC
	${T8_SRC} ${T8_CMESH_SRC} ${T8_DATA_SRC} ${T8_FOREST_SRC} ${T8_GEOMETRY_SRC} ${T8_SCHEMES_SRC} ${T8_VTK_SRC}
    )
endif()

set_target_properties( T8 PROPERTIES OUTPUT_NAME t8 )

target_compile_definitions( T8 PUBLIC T8_CMAKE_BUILD )
target_compile_definitions( T8 PUBLIC T8_CC="${CMAKE_C_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CFLAGS="${CMAKE_C_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_CPP="${CMAKE_CXX_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CPPFLAGS="${CMAKE_CXX_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LDFLAGS="${CMAKE_SHARED_LINKER_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LIBS="Not available with CMake builds" )


find_package( Git REQUIRED )
execute_process( COMMAND ${GIT_EXECUTABLE} describe --tags 
                 COMMAND cut -c 2- 
                 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                 OUTPUT_VARIABLE T8CODE_VERSION 
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
                 OUTPUT_VARIABLE T8CODE_VERSION_POINT 
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

target_compile_definitions( T8 PUBLIC T8_PACKAGE_STRING="t8 ${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION="${T8CODE_VERSION}" )
target_compile_definitions( T8 PUBLIC T8_VERSION_MAJOR=${T8CODE_VERSION_MAJOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_MINOR=${T8CODE_VERSION_MINOR} )
target_compile_definitions( T8 PUBLIC T8_VERSION_POINT=${T8CODE_VERSION_POINT} )


target_include_directories( T8 PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_link_libraries( T8 PUBLIC P4EST::P4EST SC::SC )

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_DEBUG )
endif()

if ( T8CODE_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPIIO )
    target_link_libraries( T8 PUBLIC MPI::MPI_C )
endif()

if( T8CODE_ENABLE_VTK )
    target_compile_definitions( T8 PUBLIC T8_VTK_VERSION_USED="${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" )
    target_compile_definitions( T8 PUBLIC T8_WITH_VTK=1 )
    target_include_directories( T8 PUBLIC ${VTK_INCLUDE_DIRS} )
    target_link_libraries( T8 PUBLIC ${VTK_LIBRARIES} )
endif()

install( FILES ${T8_PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include )
install( TARGETS T8 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

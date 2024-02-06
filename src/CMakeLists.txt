add_library( T8 )
target_include_directories( T8 PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_link_libraries( T8 PUBLIC P4EST::P4EST SC::SC $<$<BOOL:${T8CODE_ENABLE_MPI}>:MPI::MPI_C> )
target_sources( T8 PRIVATE 
    t8.c 
    t8_eclass.c 
    t8_mesh.c 
    t8_element.c 
    t8_element_cxx.cxx 
    t8_element_c_interface.cxx 
    t8_refcount.c 
    t8_cmesh/t8_cmesh.c 
    t8_cmesh/t8_cmesh_occ.cxx 
    t8_cmesh/t8_cmesh_cxx.cxx 
    t8_cmesh/t8_cmesh_triangle.c 
    t8_cmesh/t8_cmesh_vtk_writer.c 
    t8_cmesh/t8_cmesh_stash.c 
    t8_cmesh/t8_cmesh_vtk_reader.cxx 
    t8_cmesh/t8_cmesh_save.c 
    t8_cmesh/t8_cmesh_netcdf.c 
    t8_cmesh/t8_cmesh_trees.c 
    t8_cmesh/t8_cmesh_commit.c 
    t8_cmesh/t8_cmesh_partition.c 
    t8_cmesh/t8_cmesh_copy.c 
    t8_data/t8_shmem.c 
    t8_cmesh/t8_cmesh_geometry.cxx 
    t8_cmesh/t8_cmesh_examples.c 
    t8_cmesh/t8_cmesh_helpers.c 
    t8_data/t8_containers.cxx 
    t8_cmesh/t8_cmesh_offset.c 
    t8_cmesh/t8_cmesh_readmshfile.cxx 
    t8_forest/t8_forest.c 
    t8_forest/t8_forest_adapt.cxx 
    t8_geometry/t8_geometry.cxx 
    t8_geometry/t8_geometry_helpers.c 
    t8_geometry/t8_geometry_base.cxx 
    t8_geometry/t8_geometry_with_vertices.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_analytic.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_occ.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_linear.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_zero.cxx 
    t8_geometry/t8_geometry_implementations/t8_geometry_examples.cxx 
    t8_forest/t8_forest_partition.cxx 
    t8_forest/t8_forest_cxx.cxx 
    t8_forest/t8_forest_private.c 
    t8_forest/t8_forest_vtk.cxx 
    t8_forest/t8_forest_ghost.cxx 
    t8_forest/t8_forest_iterate.cxx 
    t8_version.c 
    t8_vtk.c 
    t8_forest/t8_forest_balance.cxx 
    t8_forest/t8_forest_netcdf.cxx 
    t8_element_shape.c 
    t8_netcdf.c 
    t8_cmesh/t8_cmesh_testcases.c 
    t8_vtk/t8_vtk_polydata.cxx 
    t8_vtk/t8_vtk_unstructured.cxx 
    t8_vtk/t8_vtk_parallel.cxx 
    t8_vtk/t8_vtk_reader.cxx 
    t8_schemes/t8_default/t8_default_cxx.cxx
    t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.cxx
    t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.cxx
    t8_schemes/t8_default/t8_default_hex/t8_dhex_bits.c
    t8_schemes/t8_default/t8_default_line/t8_default_line_cxx.cxx
    t8_schemes/t8_default/t8_default_line/t8_dline_bits.c
    t8_schemes/t8_default/t8_default_prism/t8_default_prism_cxx.cxx
    t8_schemes/t8_default/t8_default_prism/t8_dprism_bits.c
    t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.cxx
    t8_schemes/t8_default/t8_default_quad/t8_dquad_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_default_tet_cxx.cxx
    t8_schemes/t8_default/t8_default_tet/t8_dtet_bits.c
    t8_schemes/t8_default/t8_default_tet/t8_dtet_connectivity.c
    t8_schemes/t8_default/t8_default_tri/t8_default_tri_cxx.cxx
    t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.c
    t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.c
    t8_schemes/t8_default/t8_default_vertex/t8_default_vertex_cxx.cxx
    t8_schemes/t8_default/t8_default_vertex/t8_dvertex_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_default_pyramid_cxx.cxx
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.c
    t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_connectivity.c
)

target_compile_definitions( T8 PUBLIC T8_CMAKE_BUILD )
target_compile_definitions( T8 PUBLIC T8_CC="${CMAKE_C_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CFLAGS="${CMAKE_C_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_CPP="${CMAKE_CXX_COMPILER}" )
target_compile_definitions( T8 PUBLIC T8_CPPFLAGS="${CMAKE_CXX_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LDFLAGS="${CMAKE_SHARED_LINKER_FLAGS}" )
target_compile_definitions( T8 PUBLIC T8_LIBS="???" )

target_compile_definitions( T8 PUBLIC T8_PACKAGE_STRING="t8 1.6.1.???" )
target_compile_definitions( T8 PUBLIC T8_VERSION="1.6.1.???" )
target_compile_definitions( T8 PUBLIC T8_VERSION_MAJOR=1 )
target_compile_definitions( T8 PUBLIC T8_VERSION_MINOR=6 )
target_compile_definitions( T8 PUBLIC T8_VERSION_POINT=1.??? )

if ( T8CODE_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPI )
    target_compile_definitions( T8 PUBLIC T8_ENABLE_MPIIO )
endif()
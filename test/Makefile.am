#This file is part of t8code
#Non-recursive Makefile.am in test
#Included from toplevel directory

t8code_googletestincludedir=

t8code_googletest_installed_headers = \
  test/t8_gtest_custom_assertion.hxx

t8code_googletest_internal_headers = \
  thirdparty/googletest-mpi/gtest/gtest.h \
  test/t8_gtest_macros.hxx \
  test/t8_schemes/t8_gtest_dfs_base.hxx \
  test/t8_cmesh_generator/t8_cmesh_example_sets.hxx \
  test/t8_cmesh_generator/t8_gtest_cmesh_cartestian_product.hxx \
  test/t8_cmesh_generator/t8_gtest_cmesh_sum_of_sets.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_bigmesh_param.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_comm.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_disjoint_bricks_param.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_empty.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_from_class_param.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_hypercube_pad.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_hypercube_param.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_prism_cake_param.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_params.hxx \
  test/t8_cmesh_generator/t8_cmesh_parametrized_examples/t8_cmesh_new_periodic.hxx
  
t8code_googletest_programs = \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_lagrange \
  test/t8_gtest_cmesh_bcast \
  test/t8_schemes/t8_gtest_nca \
  test/t8_schemes/t8_gtest_pyra_connectivity \
  test/t8_schemes/t8_gtest_face_neigh \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_linear \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_cad \
  test/t8_gtest_eclass \
  test/t8_gtest_vec \
  test/t8_gtest_mat \
  test/t8_gtest_refcount \
  test/t8_gtest_occ_linkage \
  test/t8_gtest_version \
  test/t8_schemes/t8_gtest_init_linear_id \
  test/t8_gtest_basics \
  test/t8_schemes/t8_gtest_ancestor \
  test/t8_cmesh/t8_gtest_hypercube \
  test/t8_schemes/t8_gtest_element_count_leaves \
  test/t8_schemes/t8_gtest_element_ref_coords \
  test/t8_geometry/t8_gtest_geometry_handling \
  test/t8_schemes/t8_gtest_descendant \
  test/t8_schemes/t8_gtest_find_parent \
  test/t8_schemes/t8_gtest_equal \
  test/t8_schemes/t8_gtest_root \
  test/t8_cmesh/t8_gtest_cmesh_face_is_boundary \
  test/t8_cmesh/t8_gtest_cmesh_partition \
  test/t8_cmesh/t8_gtest_cmesh_copy \
  test/t8_cmesh/t8_gtest_cmesh_set_partition_offsets \
  test/t8_cmesh/t8_gtest_cmesh_set_join_by_vertices \
  test/t8_forest/t8_gtest_element_volume \
  test/t8_cmesh/t8_gtest_multiple_attributes \
  test/t8_cmesh/t8_gtest_cmesh_add_attributes_when_derive \
  test/t8_cmesh/t8_gtest_attribute_gloidx_array \
  test/t8_schemes/t8_gtest_successor \
  test/t8_schemes/t8_gtest_boundary_extrude \
  test/t8_forest/t8_gtest_search \
  test/t8_gtest_netcdf_linkage \
  test/t8_gtest_vtk_linkage \
  test/t8_data/t8_gtest_shmem \
  test/t8_forest/t8_gtest_half_neighbors \
  test/t8_forest/t8_gtest_find_owner \
  test/t8_forest/t8_gtest_forest_face_normal \
  test/t8_schemes/t8_gtest_face_descendant \
  test/t8_geometry/t8_gtest_point_inside \
  test/t8_forest/t8_gtest_user_data \
  test/t8_forest/t8_gtest_transform \
  test/t8_forest/t8_gtest_ghost_exchange \
  test/t8_forest/t8_gtest_ghost_delete \
  test/t8_forest/t8_gtest_ghost_and_owner \
  test/t8_forest/t8_gtest_forest_commit \
  test/t8_forest/t8_gtest_balance \
  test/t8_IO/t8_gtest_vtk_reader \
  test/t8_forest_incomplete/t8_gtest_permute_hole \
  test/t8_forest_incomplete/t8_gtest_recursive \
  test/t8_forest_incomplete/t8_gtest_iterate_replace \
  test/t8_forest_incomplete/t8_gtest_empty_local_tree \
  test/t8_forest_incomplete/t8_gtest_empty_global_tree \
  test/t8_cmesh/t8_gtest_cmesh_tree_vertices_negative_volume \
  test/t8_schemes/t8_gtest_default \
  test/t8_schemes/t8_gtest_pack_unpack \
  test/t8_schemes/t8_gtest_child_parent_face \
  test/t8_cmesh_generator/t8_gtest_cmesh_generator_test \
  test/t8_forest/t8_gtest_build_forest_from_geo_domain

test_t8_IO_t8_gtest_vtk_reader_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_IO/t8_gtest_vtk_reader.cxx

test_t8_gtest_cmesh_bcast_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_bcast.cxx

test_t8_schemes_t8_gtest_nca_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_nca.cxx

test_t8_schemes_t8_gtest_pyra_connectivity_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_pyra_connectivity.cxx

test_t8_schemes_t8_gtest_face_neigh_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_face_neigh.cxx

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_linear_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_linear.cxx

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_cad_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_cad.cxx

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_lagrange_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_geometry/t8_geometry_implementations/t8_gtest_geometry_lagrange.cxx

test_t8_gtest_eclass_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_eclass.cxx

test_t8_gtest_vec_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_vec.cxx

test_t8_gtest_mat_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_mat.cxx

test_t8_gtest_refcount_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_refcount.cxx

test_t8_gtest_occ_linkage_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_occ_linkage.cxx

test_t8_gtest_version_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_version.cxx

test_t8_schemes_t8_gtest_init_linear_id_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_init_linear_id.cxx

test_t8_gtest_basics_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_basics.cxx

test_t8_schemes_t8_gtest_ancestor_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_ancestor.cxx

test_t8_cmesh_t8_gtest_hypercube_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_hypercube.cxx

test_t8_cmesh_t8_gtest_cmesh_set_join_by_vertices_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_set_join_by_vertices.cxx

test_t8_schemes_t8_gtest_element_count_leaves_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_element_count_leaves.cxx

test_t8_schemes_t8_gtest_element_ref_coords_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_element_ref_coords.cxx

test_t8_geometry_t8_gtest_geometry_handling_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_geometry/t8_gtest_geometry_handling.cxx

test_t8_schemes_t8_gtest_descendant_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_descendant.cxx

test_t8_schemes_t8_gtest_find_parent_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_find_parent.cxx

test_t8_schemes_t8_gtest_equal_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_equal.cxx

test_t8_schemes_t8_gtest_root_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_root.cxx

test_t8_cmesh_t8_gtest_cmesh_face_is_boundary_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_face_is_boundary.cxx

test_t8_cmesh_t8_gtest_cmesh_partition_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_partition.cxx

test_t8_cmesh_t8_gtest_cmesh_set_partition_offsets_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_set_partition_offsets.cxx
  
test_t8_forest_t8_gtest_element_volume_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_element_volume.cxx

test_t8_cmesh_t8_gtest_multiple_attributes_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_multiple_attributes.cxx

test_t8_cmesh_t8_gtest_cmesh_add_attributes_when_derive_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_add_attributes_when_derive.cxx

test_t8_cmesh_t8_gtest_attribute_gloidx_array_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_attribute_gloidx_array.cxx

test_t8_schemes_t8_gtest_successor_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_successor.cxx

test_t8_schemes_t8_gtest_boundary_extrude_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_boundary_extrude.cxx

test_t8_forest_t8_gtest_search_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_search.cxx

test_t8_gtest_netcdf_linkage_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_netcdf_linkage.cxx

test_t8_gtest_vtk_linkage_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_gtest_vtk_linkage.cxx

test_t8_data_t8_gtest_shmem_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_data/t8_gtest_shmem.cxx

test_t8_forest_t8_gtest_half_neighbors_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_half_neighbors.cxx

test_t8_forest_t8_gtest_find_owner_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_find_owner.cxx

test_t8_forest_t8_gtest_forest_face_normal_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_forest_face_normal.cxx

test_t8_schemes_t8_gtest_face_descendant_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_face_descendant.cxx

test_t8_geometry_t8_gtest_point_inside_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_geometry/t8_gtest_point_inside.cxx

test_t8_forest_t8_gtest_user_data_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_user_data.cxx

test_t8_forest_t8_gtest_transform_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_transform.cxx

test_t8_forest_t8_gtest_ghost_exchange_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_ghost_exchange.cxx

test_t8_forest_t8_gtest_ghost_delete_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_ghost_delete.cxx

test_t8_forest_t8_gtest_ghost_and_owner_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_ghost_and_owner.cxx

test_t8_forest_t8_gtest_forest_commit_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_forest_commit.cxx

test_t8_forest_t8_gtest_balance_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_balance.cxx

test_t8_forest_incomplete_t8_gtest_permute_hole_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest_incomplete/t8_gtest_permute_hole.cxx

test_t8_forest_incomplete_t8_gtest_recursive_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest_incomplete/t8_gtest_recursive.cxx

test_t8_forest_incomplete_t8_gtest_iterate_replace_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest_incomplete/t8_gtest_iterate_replace.cxx

test_t8_forest_incomplete_t8_gtest_empty_local_tree_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest_incomplete/t8_gtest_empty_local_tree.cxx

test_t8_forest_incomplete_t8_gtest_empty_global_tree_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest_incomplete/t8_gtest_empty_global_tree.cxx

test_t8_cmesh_t8_gtest_cmesh_tree_vertices_negative_volume_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_tree_vertices_negative_volume.cxx

test_t8_schemes_t8_gtest_default_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_default.cxx

test_t8_schemes_t8_gtest_pack_unpack_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_pack_unpack.cxx  

test_t8_schemes_t8_gtest_child_parent_face_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_schemes/t8_gtest_child_parent_face.cxx 

test_t8_cmesh_generator_t8_gtest_cmesh_generator_test_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh_generator/t8_gtest_cmesh_generator_test.cxx

test_t8_cmesh_t8_gtest_cmesh_copy_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_cmesh_copy.cxx

test_t8_forest_t8_gtest_build_forest_from_geo_domain_SOURCES = \
  test/t8_gtest_main.cxx \
  test/t8_forest/t8_gtest_build_forest_from_geo_domain.cxx

#define ld and cpp flags for all targets
t8_gtest_target_ld_add = $(LDADD) test/libgtest.la
t8_gtest_target_ld_flags = $(AM_LDFLAGS) -pthread
t8_gtest_target_cpp_flags = $(AM_CPPFLAGS) -I$(top_srcdir)/thirdparty/googletest-mpi/ \
  -pthread
#flags to add when we compile with mpi
t8_gtest_target_mpi_cpp_flags = -DGTEST_HAS_MPI=0

test_t8_gtest_cmesh_bcast_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_cmesh_bcast_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_cmesh_bcast_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_nca_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_nca_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_nca_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_pyra_connectivity_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_pyra_connectivity_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_pyra_connectivity_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_face_neigh_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_face_neigh_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_face_neigh_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_linear_LDADD = $(t8_gtest_target_ld_add)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_linear_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_linear_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_cad_LDADD = $(t8_gtest_target_ld_add)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_cad_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_cad_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_lagrange_LDADD = $(t8_gtest_target_ld_add)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_lagrange_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_lagrange_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_eclass_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_eclass_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_eclass_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_vec_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_vec_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_vec_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_mat_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_mat_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_mat_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_refcount_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_refcount_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_refcount_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_occ_linkage_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_occ_linkage_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_occ_linkage_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_version_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_version_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_version_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_init_linear_id_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_init_linear_id_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_init_linear_id_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_basics_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_basics_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_basics_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_ancestor_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_ancestor_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_ancestor_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_hypercube_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_hypercube_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_hypercube_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_set_join_by_vertices_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_set_join_by_vertices_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_set_join_by_vertices_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_element_count_leaves_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_element_count_leaves_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_element_count_leaves_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_element_ref_coords_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_element_ref_coords_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_element_ref_coords_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_geometry_t8_gtest_geometry_handling_LDADD = $(t8_gtest_target_ld_add)
test_t8_geometry_t8_gtest_geometry_handling_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_geometry_t8_gtest_geometry_handling_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_descendant_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_descendant_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_descendant_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_find_parent_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_find_parent_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_find_parent_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_equal_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_equal_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_equal_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_root_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_root_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_root_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_face_is_boundary_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_face_is_boundary_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_face_is_boundary_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_partition_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_partition_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_partition_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_set_partition_offsets_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_set_partition_offsets_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_set_partition_offsets_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_element_volume_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_element_volume_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_element_volume_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_multiple_attributes_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_multiple_attributes_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_multiple_attributes_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_add_attributes_when_derive_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_add_attributes_when_derive_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_add_attributes_when_derive_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_attribute_gloidx_array_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_attribute_gloidx_array_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_attribute_gloidx_array_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_successor_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_successor_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_successor_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_boundary_extrude_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_boundary_extrude_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_boundary_extrude_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_search_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_search_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_search_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_netcdf_linkage_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_netcdf_linkage_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_netcdf_linkage_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_gtest_vtk_linkage_LDADD = $(t8_gtest_target_ld_add)
test_t8_gtest_vtk_linkage_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_gtest_vtk_linkage_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_data_t8_gtest_shmem_LDADD = $(t8_gtest_target_ld_add)
test_t8_data_t8_gtest_shmem_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_data_t8_gtest_shmem_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_half_neighbors_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_half_neighbors_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_half_neighbors_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_find_owner_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_find_owner_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_find_owner_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_forest_face_normal_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_forest_face_normal_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_forest_face_normal_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_face_descendant_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_face_descendant_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_face_descendant_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_geometry_t8_gtest_point_inside_LDADD = $(t8_gtest_target_ld_add)
test_t8_geometry_t8_gtest_point_inside_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_geometry_t8_gtest_point_inside_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_user_data_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_user_data_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_user_data_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_transform_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_transform_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_transform_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_ghost_exchange_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_ghost_exchange_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_ghost_exchange_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_ghost_delete_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_ghost_delete_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_ghost_delete_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_ghost_and_owner_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_ghost_and_owner_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_ghost_and_owner_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_forest_commit_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_forest_commit_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_forest_commit_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_balance_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_balance_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_balance_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_IO_t8_gtest_vtk_reader_LDADD = $(t8_gtest_target_ld_add)
test_t8_IO_t8_gtest_vtk_reader_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_IO_t8_gtest_vtk_reader_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_incomplete_t8_gtest_permute_hole_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_incomplete_t8_gtest_permute_hole_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_incomplete_t8_gtest_permute_hole_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_incomplete_t8_gtest_recursive_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_incomplete_t8_gtest_recursive_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_incomplete_t8_gtest_recursive_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_incomplete_t8_gtest_iterate_replace_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_incomplete_t8_gtest_iterate_replace_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_incomplete_t8_gtest_iterate_replace_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_incomplete_t8_gtest_empty_local_tree_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_incomplete_t8_gtest_empty_local_tree_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_incomplete_t8_gtest_empty_local_tree_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_incomplete_t8_gtest_empty_global_tree_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_incomplete_t8_gtest_empty_global_tree_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_incomplete_t8_gtest_empty_global_tree_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_tree_vertices_negative_volume_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_tree_vertices_negative_volume_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_tree_vertices_negative_volume_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_default_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_default_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_default_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_pack_unpack_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_pack_unpack_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_pack_unpack_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_schemes_t8_gtest_child_parent_face_LDADD = $(t8_gtest_target_ld_add)
test_t8_schemes_t8_gtest_child_parent_face_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_schemes_t8_gtest_child_parent_face_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_generator_t8_gtest_cmesh_generator_test_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_generator_t8_gtest_cmesh_generator_test_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_generator_t8_gtest_cmesh_generator_test_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_cmesh_t8_gtest_cmesh_copy_LDADD = $(t8_gtest_target_ld_add)
test_t8_cmesh_t8_gtest_cmesh_copy_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_cmesh_t8_gtest_cmesh_copy_CPPFLAGS = $(t8_gtest_target_cpp_flags)

test_t8_forest_t8_gtest_build_forest_from_geo_domain_LDADD = $(t8_gtest_target_ld_add)
test_t8_forest_t8_gtest_build_forest_from_geo_domain_LDFLAGS = $(t8_gtest_target_ld_flags)
test_t8_forest_t8_gtest_build_forest_from_geo_domain_CPPFLAGS = $(t8_gtest_target_cpp_flags)

# If we did not configure t8code with MPI we need to build Googletest
# without MPI support.
if !T8_ENABLE_MPI
test_t8_gtest_cmesh_bcast_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_nca_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_pyra_connectivity_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_face_neigh_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_linear_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_geometry_t8_geometry_implementations_t8_gtest_geometry_cad_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_eclass_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_vec_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_mat_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_refcount_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_occ_linkage_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_version_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_init_linear_id_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_basics_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_ancestor_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_hypercube_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_set_join_by_vertices_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_element_count_leaves_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_element_ref_coords_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_geometry_t8_gtest_geometry_handling_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_descendant_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_find_parent_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_equal_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_root_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_face_is_boundary_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_partition_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_set_partition_offsets_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_element_volume_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_multiple_attributes_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_add_attributes_when_derive_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_attribute_gloidx_array_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_successor_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_boundary_extrude_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_search_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_netcdf_linkage_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_gtest_vtk_linkage_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_data_t8_gtest_shmem_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_half_neighbors_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_find_owner_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_forest_face_normal_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_face_descendant_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_geometry_t8_gtest_point_inside_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_user_data_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_transform_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_ghost_exchange_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_ghost_delete_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_ghost_and_owner_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_forest_commit_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_balance_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_IO_t8_gtest_vtk_reader_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_incomplete_t8_gtest_permute_hole_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_incomplete_t8_gtest_recursive_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_incomplete_t8_gtest_iterate_replace_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_incomplete_t8_gtest_empty_local_tree_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_incomplete_t8_gtest_empty_global_tree_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_tree_vertices_negative_volume_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_default_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_pack_unpack_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_schemes_t8_gtest_child_parent_face_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_generator_t8_gtest_cmesh_generator_test_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_cmesh_t8_gtest_cmesh_copy_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
test_t8_forest_t8_gtest_build_forest_from_geo_domain_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)

endif

# Build Googletest library
check_LTLIBRARIES = test/libgtest.la
test_libgtest_la_SOURCES = thirdparty/googletest-mpi/gtest/gtest-all.cc
test_libgtest_la_CPPFLAGS = -I$(top_srcdir)/thirdparty/googletest-mpi
# If we did not configure t8code with MPI we need to build Googletest
# without MPI support.
if !T8_ENABLE_MPI
test_libgtest_la_CPPFLAGS += $(t8_gtest_target_mpi_cpp_flags)
endif
test_libgtest_la_LDFLAGS = -pthread

dist_t8code_googletestinclude_HEADERS = $(t8code_googletest_installed_headers)
dist_t8code_googletestinclude_HEADERS += $(t8code_googletest_internal_headers)

TESTS += $(t8code_googletest_programs)
check_PROGRAMS += $(t8code_googletest_programs)

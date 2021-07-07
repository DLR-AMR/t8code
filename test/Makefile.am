# This file is part of t8code
# Non-recursive Makefile.am in test
# Included from toplevel directory

t8code_test_programs = \
        test/t8_test_eclass \
        test/t8_test_bcast \
        test/t8_test_hypercube \
        test/t8_test_cmesh_copy \
        test/t8_test_cmesh_partition \
        test/t8_test_find_owner \
        test/t8_test_ghost_exchange \
        test/t8_test_ghost_and_owner \
	test/t8_test_forest_commit \
	test/t8_test_transform \
	test/t8_test_half_neighbors \
	test/t8_test_point_inside \
	test/t8_test_element_count_leafs \
	test/t8_test_search \
	test/t8_test_find_parent \
        test/t8_test_successor \
        test/t8_test_init_linear_id \
        test/t8_test_descendant \
	test/t8_test_cmesh_face_is_boundary \
	test/t8_test_element_general_function \
	test/t8_test_cmesh_readmshfile \
	test/t8_test_netcdf_linkage \
	test/t8_test_vtk_linkage \
	test/t8_test_user_data  \
        test/t8_test_pyra_face_neigh \
        test/t8_test_pyra_face_descendant

test_t8_test_eclass_SOURCES = test/t8_test_eclass.c
test_t8_test_bcast_SOURCES = test/t8_test_bcast.c
test_t8_test_hypercube_SOURCES = test/t8_test_hypercube.c
test_t8_test_cmesh_copy_SOURCES = test/t8_test_cmesh_copy.c
test_t8_test_cmesh_partition_SOURCES = test/t8_test_cmesh_partition.cxx
test_t8_test_find_owner_SOURCES = test/t8_test_find_owner.cxx
test_t8_test_ghost_exchange_SOURCES = test/t8_test_ghost_exchange.cxx
test_t8_test_ghost_and_owner_SOURCES = test/t8_test_ghost_and_owner.cxx
test_t8_test_forest_commit_SOURCES = test/t8_test_forest_commit.cxx
test_t8_test_transform_SOURCES = test/t8_test_transform.cxx
test_t8_test_half_neighbors_SOURCES = test/t8_test_half_neighbors.cxx
test_t8_test_element_count_leafs_SOURCES = test/t8_test_element_count_leafs.cxx
test_t8_test_search_SOURCES = test/t8_test_search.cxx
test_t8_test_cmesh_face_is_boundary_SOURCES = test/t8_test_cmesh_face_is_boundary.cxx
test_t8_test_element_general_function_SOURCES = test/t8_test_element_general_function.cxx
test_t8_test_point_inside_SOURCES = test/t8_test_point_inside.cxx
test_t8_test_find_parent_SOURCES = test/t8_test_find_parent.cxx
test_t8_test_successor_SOURCES = test/t8_test_successor.cxx
test_t8_test_init_linear_id_SOURCES = test/t8_test_init_linear_id.cxx
test_t8_test_descendant_SOURCES = test/t8_test_descendant.cxx
test_t8_test_cmesh_face_is_boundary_SOURCES = test/t8_test_cmesh_face_is_boundary.cxx
test_t8_test_element_general_function_SOURCES = test/t8_test_element_general_function.cxx
test_t8_test_cmesh_readmshfile_SOURCES = test/t8_test_cmesh_readmshfile.c
test_t8_test_netcdf_linkage_SOURCES = test/t8_test_netcdf_linkage.c
test_t8_test_vtk_linkage_SOURCES = test/t8_test_vtk_linkage.cxx
test_t8_test_user_data_SOURCES = test/t8_test_user_data.cxx
test_t8_test_pyra_face_neigh_SOURCES = test/t8_test_pyra_face_neigh.cxx
test_t8_test_pyra_face_descendant_SOURCES = test/t8_test_pyra_face_descendant.cxx

TESTS += $(t8code_test_programs)
check_PROGRAMS += $(t8code_test_programs)

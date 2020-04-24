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
	test/t8_test_element_count_leafs \
	test/t8_test_search \
	test/t8_test_find_parent \
	test/t8_test_cmesh_face_is_boundary

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
test_t8_test_find_parent_SOURCES = test/t8_test_find_parent.cpp
test_t8_test_cmesh_face_is_boundary_SOURCES = test/t8_test_cmesh_face_is_boundary.cxx

TESTS += $(t8code_test_programs)
check_PROGRAMS += $(t8code_test_programs)

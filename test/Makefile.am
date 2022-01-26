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
	test/t8_test_geometry \
        test/t8_test_ghost_exchange \
        test/t8_test_ghost_and_owner \
	test/t8_test_forest_commit \
	test/t8_test_transform \
	test/t8_test_half_neighbors \
	test/t8_test_point_inside \
	test/t8_test_element_count_leafs \
	test/t8_test_search \
	test/t8_test_find_parent \
	test/t8_test_cmesh_face_is_boundary \
	test/t8_test_element_general_function \
	test/t8_test_cmesh_readmshfile \
	test/t8_test_netcdf_linkage \
	test/t8_test_vtk_linkage \
	test/t8_test_shmem \
	test/t8_test_user_data

test_t8_test_eclass_SOURCES = test/t8_test_eclass.c
test_t8_test_bcast_SOURCES = test/t8_test_bcast.c
test_t8_test_hypercube_SOURCES = test/t8_test_hypercube.c
test_t8_test_cmesh_copy_SOURCES = test/t8_test_cmesh_copy.c
test_t8_test_cmesh_partition_SOURCES = test/t8_test_cmesh_partition.cxx
test_t8_test_find_owner_SOURCES = test/t8_test_find_owner.cxx
test_t8_test_geometry_SOURCES = test/t8_test_geometry.cxx
test_t8_test_ghost_exchange_SOURCES = test/t8_test_ghost_exchange.cxx
test_t8_test_ghost_and_owner_SOURCES = test/t8_test_ghost_and_owner.cxx
test_t8_test_forest_commit_SOURCES = test/t8_test_forest_commit.cxx
test_t8_test_transform_SOURCES = test/t8_test_transform.cxx
test_t8_test_half_neighbors_SOURCES = test/t8_test_half_neighbors.cxx
test_t8_test_element_count_leafs_SOURCES = test/t8_test_element_count_leafs.cxx
test_t8_test_search_SOURCES = test/t8_test_search.cxx
test_t8_test_point_inside_SOURCES = test/t8_test_point_inside.cxx
test_t8_test_find_parent_SOURCES = test/t8_test_find_parent.cxx
test_t8_test_cmesh_face_is_boundary_SOURCES = test/t8_test_cmesh_face_is_boundary.cxx
test_t8_test_element_general_function_SOURCES = test/t8_test_element_general_function.cxx
test_t8_test_cmesh_readmshfile_SOURCES = test/t8_test_cmesh_readmshfile.c
test_t8_test_netcdf_linkage_SOURCES = test/t8_test_netcdf_linkage.c
test_t8_test_vtk_linkage_SOURCES = test/t8_test_vtk_linkage.cxx
test_t8_test_shmem_SOURCES = test/t8_test_shmem.cxx
test_t8_test_user_data_SOURCES = test/t8_test_user_data.cxx

t8code_googletest_programs = \
  test/t8_gtest_main

test_t8_gtest_main_SOURCES = test/t8_gtest_main.cxx \
  test/t8_gtest_eclass.cxx \
  test/t8_gtest_vec.cxx \
  test/t8_gtest_refcount.cxx
test_t8_gtest_main_LDADD = $(LDADD) test/libgtest.la
test_t8_gtest_main_LDFLAGS = $(AM_LDFLAGS) -pthread
test_t8_gtest_main_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/thirdparty/googletest-mpi/ \
	-pthread \
  -std=c++11
# If we did not configure t8code with MPI we need to build Googletest
# without MPI support.
if !T8_ENABLE_MPI
test_t8_gtest_main_CPPFLAGS += -DGTEST_HAS_MPI=0
endif

# Build Googletest library
check_LTLIBRARIES = test/libgtest.la
test_libgtest_la_SOURCES = thirdparty/googletest-mpi/gtest/gtest-all.cc
test_libgtest_la_CPPFLAGS = -I$(top_srcdir)/thirdparty/googletest-mpi -std=c++11
# If we did not configure t8code with MPI we need to build Googletest
# without MPI support.
if !T8_ENABLE_MPI
test_libgtest_la_CPPFLAGS += -DGTEST_HAS_MPI=0
endif
test_libgtest_la_LDFLAGS = -pthread


TESTS += $(t8code_test_programs) $(t8code_googletest_programs)
check_PROGRAMS += $(t8code_test_programs) $(t8code_googletest_programs)

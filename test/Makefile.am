# This file is part of t8code
# Non-recursive Makefile.am in test
# Included from toplevel directory

t8code_test_programs = \
    test/t8_cmesh/t8_test_cmesh_partition \
    test/t8_cmesh/t8_test_cmesh_face_is_boundary \
    test/t8_cmesh/t8_test_cmesh_readmshfile \
    test/t8_forest/t8_test_find_owner \
    test/t8_forest/t8_test_ghost_exchange \
    test/t8_forest/t8_test_ghost_and_owner \
    test/t8_forest/t8_test_forest_commit \
    test/t8_forest/t8_test_transform \
    test/t8_forest/t8_test_half_neighbors \
    test/t8_forest/t8_test_search \
    test/t8_forest/t8_test_element_general_function \
    test/t8_forest/t8_test_user_data  \
    test/t8_geometry/t8_test_geometry \
    test/t8_geometry/t8_test_point_inside \
    test/t8_schemes/t8_test_element_count_leafs \
    test/t8_schemes/t8_test_find_parent \
    test/t8_schemes/t8_test_successor \
    test/t8_schemes/t8_test_descendant \
    test/t8_schemes/t8_test_pyra_face_neigh \
    test/t8_schemes/t8_test_pyra_face_descendant \
    test/t8_data/t8_test_shmem \
    test/t8_test_netcdf_linkage \
    test/t8_test_vtk_linkage \
    test/t8_forest_incomplete/t8_test_elements_remove \
	  test/t8_forest_incomplete/t8_test_remove_rand \
	  test/t8_forest_incomplete/t8_test_elements_interpolation \
	  test/t8_forest_incomplete/t8_test_recursive

test_t8_cmesh_t8_test_cmesh_copy_SOURCES = test/t8_cmesh/t8_test_cmesh_copy.c
test_t8_cmesh_t8_test_cmesh_partition_SOURCES = test/t8_cmesh/t8_test_cmesh_partition.cxx
test_t8_cmesh_t8_test_cmesh_face_is_boundary_SOURCES = test/t8_cmesh/t8_test_cmesh_face_is_boundary.cxx
test_t8_cmesh_t8_test_cmesh_readmshfile_SOURCES = test/t8_cmesh/t8_test_cmesh_readmshfile.cxx

test_t8_forest_t8_test_find_owner_SOURCES = test/t8_forest/t8_test_find_owner.cxx
test_t8_forest_t8_test_ghost_exchange_SOURCES = test/t8_forest/t8_test_ghost_exchange.cxx
test_t8_forest_t8_test_ghost_and_owner_SOURCES = test/t8_forest/t8_test_ghost_and_owner.cxx
test_t8_forest_t8_test_forest_commit_SOURCES = test/t8_forest/t8_test_forest_commit.cxx
test_t8_forest_t8_test_transform_SOURCES = test/t8_forest/t8_test_transform.cxx
test_t8_forest_t8_test_half_neighbors_SOURCES = test/t8_forest/t8_test_half_neighbors.cxx
test_t8_forest_t8_test_search_SOURCES = test/t8_forest/t8_test_search.cxx
test_t8_forest_t8_test_element_general_function_SOURCES = test/t8_forest/t8_test_element_general_function.cxx
test_t8_forest_t8_test_user_data_SOURCES = test/t8_forest/t8_test_user_data.cxx

test_t8_schemes_t8_test_pyra_face_neigh_SOURCES = test/t8_schemes/t8_test_pyra_face_neigh.cxx
test_t8_schemes_t8_test_pyra_face_descendant_SOURCES = test/t8_schemes/t8_test_pyra_face_descendant.cxx
test_t8_schemes_t8_test_find_parent_SOURCES = test/t8_schemes/t8_test_find_parent.cxx
test_t8_schemes_t8_test_successor_SOURCES = test/t8_schemes/t8_test_successor.cxx
test_t8_schemes_t8_test_descendant_SOURCES = test/t8_schemes/t8_test_descendant.cxx
test_t8_schemes_t8_test_element_count_leafs_SOURCES = test/t8_schemes/t8_test_element_count_leafs.cxx

test_t8_geometry_t8_test_point_inside_SOURCES = test/t8_geometry/t8_test_point_inside.cxx
test_t8_geometry_t8_test_geometry_SOURCES = test/t8_geometry/t8_test_geometry.cxx

test_t8_data_t8_test_shmem_SOURCES = test/t8_data/t8_test_shmem.cxx

test_t8_test_netcdf_linkage_SOURCES = test/t8_test_netcdf_linkage.c
test_t8_test_vtk_linkage_SOURCES = test/t8_test_vtk_linkage.cxx

test_t8_forest_incomplete_t8_test_elements_remove_SOURCES = test/t8_forest_incomplete/t8_test_elements_remove.cxx
test_t8_forest_incomplete_t8_test_remove_rand_SOURCES = test/t8_forest_incomplete/t8_test_remove_rand.cxx
test_t8_forest_incomplete_t8_test_elements_interpolation_SOURCES = test/t8_forest_incomplete/t8_test_elements_interpolation.cxx
test_t8_forest_incomplete_t8_test_recursive_SOURCES = test/t8_forest_incomplete/t8_test_recursive.cxx

t8code_googletest_programs = \
  test/t8_gtest_main

test_t8_gtest_main_SOURCES = test/t8_gtest_main.cxx \
  test/t8_cmesh/t8_gtest_bcast.cxx \
  test/t8_schemes/t8_gtest_nca.cxx \
  test/t8_schemes/t8_gtest_pyra_connectivity.cxx \
  test/t8_geometry/t8_gtest_geometry_occ.cxx \
  test/t8_gtest_eclass.cxx \
  test/t8_gtest_vec.cxx \
  test/t8_gtest_refcount.cxx \
  test/t8_gtest_occ_linkage.cxx \
  test/t8_gtest_version.cxx \
  test/t8_schemes/t8_gtest_init_linear_id.cxx \
  test/t8_gtest_basics.cxx \
  test/t8_schemes/t8_gtest_ancestor.cxx \
  test/t8_cmesh/t8_gtest_hypercube.cxx \
  test/t8_cmesh/t8_gtest_cmesh_copy.cxx

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

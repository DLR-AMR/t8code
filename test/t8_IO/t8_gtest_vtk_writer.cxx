/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2024 the developers

t8code is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

t8code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with t8code; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <gtest/gtest.h>
#include <t8_vtk/t8_vtk_writer.hxx>
#include <t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_vtk/t8_vtk_writer_c_interface.h>

/**
 * Create a hybrid forest or a cmesh
 * 
 * \tparam grid_t 
 * \return grid_t 
 */
template <typename grid_t>
grid_t
make_grid ();

template <>
t8_cmesh_t
make_grid<t8_cmesh_t> ()
{
  return t8_cmesh_new_full_hybrid (sc_MPI_COMM_WORLD);
}
template <>
t8_forest_t
make_grid<t8_forest_t> ()
{
  t8_cmesh_t cmesh = make_grid<t8_cmesh_t> ();
  t8_scheme_cxx *scheme = t8_scheme_new_default_cxx ();
  return t8_forest_new_uniform (cmesh, scheme, 2, 0, sc_MPI_COMM_WORLD);
}

/**
 * Destroy a forest or a cmesh
 * 
 * \tparam grid_t 
 * \param[in] grid A pointer to a forest/cmesh that will get unrefed. 
 */
template <typename grid_t>
void
destroy_grid (grid_t *grid);

template <>
void
destroy_grid<t8_cmesh_t> (t8_cmesh_t *cmesh)
{
  t8_cmesh_unref (cmesh);
}

template <>
void
destroy_grid<t8_forest_t> (t8_forest_t *forest)
{
  t8_forest_unref (forest);
}

template <typename grid_t>
int
use_c_interface (const grid_t grid, const char *fileprefix, const int write_treeid, const int write_mpirank,
                 const int write_level, const int write_element_id, const int curved_flag, const int write_ghosts,
                 const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm);

template <>
int
use_c_interface<t8_forest_t> (const t8_forest_t grid, const char *fileprefix, const int write_treeid,
                              const int write_mpirank, const int write_level, const int write_element_id,
                              const int curved_flag, const int write_ghosts, const int num_data,
                              t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  return t8_forest_vtk_write_file_via_API (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                           curved_flag, write_ghosts, num_data, data);
#else
  return t8_forest_vtk_write_file (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                   curved_flag, write_ghosts, num_data, data);
#endif
}

template <>
int
use_c_interface<t8_cmesh_t> (const t8_cmesh_t grid, const char *fileprefix, const int write_treeid,
                             const int write_mpirank, const int write_level, const int write_element_id,
                             const int curved_flag, const int write_ghosts, const int num_data,
                             t8_vtk_data_field_t *data, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  return t8_cmesh_vtk_write_file_via_API (grid, fileprefix, comm);
#else
  return t8_cmesh_vtk_write_file (grid, fileprefix);
#endif
}

/**
 * Templated class to test the vtk writer for forests and cmeshes. 
 * 
 * @tparam grid_t 
 */
template <typename grid_t>
class vtk_writer_test: public testing::Test {
 protected:
  void
  SetUp () override
  {
    grid = make_grid<grid_t> ();
    writer = new vtk_writer<grid_t> (true, true, true, true, true, true, std::string ("test_vtk"), 0, NULL,
                                     sc_MPI_COMM_WORLD);
  }

  int
  grid_c_interface ()
  {
    return use_c_interface (grid, "test_vtk_c_interface", 1, 1, 1, 1, 1, 1, 0, NULL, sc_MPI_COMM_WORLD);
  }

  void
  TearDown () override
  {
    destroy_grid (&grid);
  }

  grid_t grid;
  vtk_writer<grid_t> *writer;
};

TYPED_TEST_SUITE_P (vtk_writer_test);

/**
 * Test the API-writer and check if the output was created successfully. 
 * 
 */
TYPED_TEST_P (vtk_writer_test, write_vtk)
{
#if T8_WITH_VTK
  EXPECT_TRUE (this->writer->write_with_API (this->grid));
#else
  EXPECT_TRUE (this->writer->write_ASCII (this->grid));
#endif
}

TYPED_TEST_P (vtk_writer_test, c_interface)
{
  EXPECT_TRUE (this->grid_c_interface ());
}

REGISTER_TYPED_TEST_SUITE_P (vtk_writer_test, write_vtk, c_interface);

using GridTypes = ::testing::Types<t8_cmesh_t, t8_forest_t>;

INSTANTIATE_TYPED_TEST_SUITE_P (Test_vtk_writer, vtk_writer_test, GridTypes, );
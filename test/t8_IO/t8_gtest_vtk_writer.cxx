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
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_vtk/t8_vtk_writer.h>

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
  const t8_scheme *scheme = t8_scheme_new_default ();
  const int uniform_level = 2;
  const int do_ghosts = 1;
  return t8_forest_new_uniform (cmesh, scheme, uniform_level, do_ghosts, sc_MPI_COMM_WORLD);
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
static int
use_c_interface (const grid_t grid, const char *fileprefix, const int write_treeid, const int write_mpirank,
                 const int write_level, const int write_element_id, const int curved_flag, const int write_ghosts,
                 const int num_data, t8_vtk_data_field_t *data, sc_MPI_Comm comm, const bool use_api);

template <>
int
use_c_interface<t8_forest_t> (const t8_forest_t grid, const char *fileprefix, const int write_treeid,
                              const int write_mpirank, const int write_level, const int write_element_id,
                              [[maybe_unused]] const int curved_flag, const int write_ghosts, const int num_data,
                              t8_vtk_data_field_t *data, [[maybe_unused]] sc_MPI_Comm comm, const bool use_api)
{
  if (use_api)
  {
    #if T8_ENABLE_VTK
      return t8_forest_vtk_write_file_via_API (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                           curved_flag, write_ghosts, num_data, data);
    #endif
      t8_errorf ("Testing API without linking against VTK, falling back to non-VTK API\n");
  }
    return t8_forest_vtk_write_file (grid, fileprefix, write_treeid, write_mpirank, write_level, write_element_id,
                                   write_ghosts, num_data, data);
}

template <>
int
use_c_interface<t8_cmesh_t> (const t8_cmesh_t grid, const char *fileprefix, [[maybe_unused]] const int write_treeid,
                             [[maybe_unused]] const int write_mpirank, [[maybe_unused]] const int write_level,
                             [[maybe_unused]] const int write_element_id, [[maybe_unused]] const int curved_flag,
                             [[maybe_unused]] const int write_ghosts, [[maybe_unused]] const int num_data,
                             [[maybe_unused]] t8_vtk_data_field_t *data, [[maybe_unused]] sc_MPI_Comm comm, const bool use_api)
{
  if (use_api)
  {
    #if T8_ENABLE_VTK
      return t8_cmesh_vtk_write_file_via_API (grid, fileprefix, comm);
    #endif
    t8_errorf ("Testing API without linking against VTK, falling back to non-VTK API\n");
  }
  return t8_cmesh_vtk_write_file (grid, fileprefix);
}

void
vtk_writer_test_fill_data (const t8_locidx_t cells_to_write_count, std::vector<double> &scalar_data,
                           std::vector<double> &vector_data)
{
  scalar_data.resize (cells_to_write_count);
  vector_data.resize (3 * cells_to_write_count);
  // Fill scalar data vector with entries 0,1/10,...,(N-1)/10
  //    vector[n] = (n/10.)
  // Fill vector data vector with entries (0, 0, 42), (0.1,-0.1,42), ...
  //    vector[n] = (n/10.,-n/10., 42)
  std::generate (scalar_data.begin (), scalar_data.end (), [n = 0] () mutable { return (n++) / 10.; });
  std::generate (vector_data.begin (), vector_data.end (), [n = 0, scalar_data] () mutable {
    double scalar_value = scalar_data[n / 3];
    double vector_values[3] = { scalar_value, -scalar_value, 42. };
    return vector_values[n++ % 3];
  });
}

/**
 * Templated class to test the vtk writer for forests and cmeshes. 
 * 
 * @tparam grid_t 
 */
template <typename T>
struct vtk_writer_test: public testing::Test
{
 protected:
  void
  SetUp () override
  {
    grid = make_grid<grid_t> ();
    use_api = T::use_api;
    const t8_locidx_t cells_to_write_count = num_cells_to_write (grid, 1);

    // Fill the test data vectors with dummy data.
    vtk_writer_test_fill_data (cells_to_write_count, scalar_data, vector_data);

    // Fill the vtk_data descriptors
    vtk_data[0].type = T8_VTK_SCALAR;
    strncpy (vtk_data[0].description, "Testdata scalar i/10.", BUFSIZ);
    vtk_data[0].data = scalar_data.data ();
    vtk_data[1].type = T8_VTK_VECTOR;
    strncpy (vtk_data[1].description, "Testdata vector (i/10.,-i/10.,42)", BUFSIZ);
    vtk_data[1].data = vector_data.data ();

    writer = new vtk_writer<grid_t> (true, true, true, true, true, true, std::string ("test_vtk_writer"), num_vtk_data,
                                     vtk_data, sc_MPI_COMM_WORLD);
  }

  int
  grid_c_interface ()
  {
    return use_c_interface (grid, "test_vtk_writer_c_interface", 1, 1, 1, 1, 1, 1, num_vtk_data, vtk_data,
                            sc_MPI_COMM_WORLD, use_api);
  }

  void
  TearDown () override
  {
    if (writer != nullptr) {
      delete writer;
      writer = nullptr;
    }
    destroy_grid (&grid);
  }
  using grid_t = typename T::grid_t;
  grid_t grid;
  vtk_writer<grid_t> *writer;
  std::vector<double> scalar_data;             // Scalar data used for data output
  std::vector<double> vector_data;             // Vector data used for data output
  static constexpr int num_vtk_data = 2;       // Number of data items
  t8_vtk_data_field_t vtk_data[num_vtk_data];  // The metadata descriptor for data output
  bool use_api = false;                        // Should we use the api or not
};

TYPED_TEST_SUITE_P (vtk_writer_test);

/**
 * Test the API-writer and check if the output was created successfully. 
 * 
 */
TYPED_TEST_P (vtk_writer_test, write_vtk)
{
  if (this->use_api)
  {
    #if T8_ENABLE_VTK
      EXPECT_TRUE (this->writer->write_with_API (this->grid));
    #endif
  }
  else
  {
    EXPECT_TRUE (this->writer->write_ASCII (this->grid));
  }
}

TYPED_TEST_P (vtk_writer_test, c_interface)
{
  EXPECT_TRUE (this->grid_c_interface ());
}

REGISTER_TYPED_TEST_SUITE_P (vtk_writer_test, write_vtk, c_interface);

template <typename T, bool api_usage>
struct TestConfig{
  using grid_t = T;
  static constexpr bool use_api = api_usage;
};
using GridTypes = ::testing::Types<
  TestConfig<t8_cmesh_t, false>,
  TestConfig<t8_cmesh_t, true>,
  TestConfig<t8_forest_t, false>,
  TestConfig<t8_forest_t, true>
>;

// Name generator for typed-tests: produces readable test-case suffixes such as
// "t8_cmesh_noapi", "t8_forest_api", ...
struct GridTypesNameGenerator {
  template <typename T>
  static std::string GetName(int)
  {
    // use std::is_same without if constexpr for broader compiler compatibility
    std::string name = std::is_same<typename T::grid_t, t8_cmesh_t>::value ? "t8_cmesh" : "t8_forest";
    name += T::use_api ? "_with_api" : "_no_api";
    return name;
  }
};

INSTANTIATE_TYPED_TEST_SUITE_P (Test_vtk_writer, vtk_writer_test, GridTypes, GridTypesNameGenerator);

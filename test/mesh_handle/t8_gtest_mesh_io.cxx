/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/**
 * \file t8_gtest_mesh_io.cxx
 * Tests the vtk output functions of the mesh handle defined in mesh_io.hxx.
 * Besides checking that the output is successful, the output of the handle is compared to the output of the forest.
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8.h>

#include <mesh_handle/mesh_handle.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_io.h>
#include <t8_vtk/t8_vtk.h>
#include <t8_vector/t8_vec.hxx>
#include <cstring>
#include <fstream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include <cstdio>
#include <filesystem>

/** Parametrized test fixture for the mesh io tests.
 * Provides a mesh with ghosts and exemplary per element data for all local and ghost elements.
 */
struct t8_mesh_io_test: public testing::TestWithParam<std::tuple<t8_eclass_t, int>>
{
 protected:
  using mesh_class = t8_mesh_handle::mesh<>; /**< The mesh handle class used for testing. */

  void
  SetUp () override
  {
    // All vtk files of a test are written into output_dir, which is deleted again in TearDown.
    int mpiret = sc_MPI_Comm_rank (sc_MPI_COMM_WORLD, &mpirank);
    SC_CHECK_MPI (mpiret);
    if (mpirank == 0) {
      std::filesystem::create_directories (output_dir);
    }
    //All processes have to wait until the directory exists before they write into it.
    mpiret = sc_MPI_Barrier (sc_MPI_COMM_WORLD);
    SC_CHECK_MPI (mpiret);
    const t8_eclass_t eclass = std::get<0> (GetParam ());
    const int level = std::get<1> (GetParam ());
    // Create the mesh with ghosts to be able to test the output of ghost elements.
    mesh = t8_mesh_handle::handle_hypercube_uniform_default<mesh_class> (eclass, level, sc_MPI_COMM_WORLD, true, true,
                                                                         false);
    num_local_elements = mesh->get_num_local_elements ();
    num_ghosts = mesh->get_num_ghosts ();

    // Use the level, the volume and the centroid of each local and ghost element as data to write.
    for (t8_locidx_t ielem = 0; ielem < num_local_elements + num_ghosts; ++ielem) {
      levels.push_back ((*mesh)[ielem].get_level ());
      volumes.push_back ((*mesh)[ielem].get_volume ());
      const t8_3D_vec centroid = (*mesh)[ielem].get_centroid ();
      centroids.push_back (centroid);
      centroid_components.insert (centroid_components.end (), centroid.begin (), centroid.end ());
    }

    // Fill the data field descriptors as expected by the forest vtk writer. All scalar fields have to come first.
    vtk_data[0].type = T8_VTK_SCALAR;
    strncpy (vtk_data[0].description, "level", T8_VTK_MAX_STRING_LENGTH);
    vtk_data[0].data = levels.data ();
    vtk_data[1].type = T8_VTK_SCALAR;
    strncpy (vtk_data[1].description, "volume", T8_VTK_MAX_STRING_LENGTH);
    vtk_data[1].data = volumes.data ();
    vtk_data[2].type = T8_VTK_VECTOR;
    strncpy (vtk_data[2].description, "centroid", T8_VTK_MAX_STRING_LENGTH);
    vtk_data[2].data = centroid_components.data ();
  }

  void
  TearDown () override
  {
    // All processes have to be done writing and reading before the files are deleted.
    const int mpiret = sc_MPI_Barrier (sc_MPI_COMM_WORLD);
    SC_CHECK_MPI (mpiret);
    if (mpirank == 0) {
      std::filesystem::remove_all (output_dir);
    }
  }

  /** Get the file prefix for the vtk output \a name inside \a output_dir. */
  std::string
  output_prefix (const std::string &name) const
  {
    return output_dir + "/" + name;
  }

  const std::string output_dir = "t8_gtest_mesh_io_output"; /**< Directory for all files written. */
  int mpirank;                                              /**< The rank of this process. */
  std::unique_ptr<mesh_class> mesh;                         /**< The mesh to write. */
  t8_locidx_t num_local_elements;                           /**< Number of local elements of \a mesh. */
  t8_locidx_t num_ghosts;                                   /**< Number of ghost elements of \a mesh. */
  std::vector<double> levels;                               /**< Level of each local and ghost element. */
  std::vector<double> volumes;                              /**< Volume of each local and ghost element. */
  std::vector<t8_3D_vec> centroids;                         /**< Centroid of each local and ghost element. */
  std::vector<double> centroid_components;    /**< The centroids as flat array as expected by the forest. */
  static constexpr int num_vtk_data = 3;      /**< Number of data fields in \a vtk_data. */
  t8_vtk_data_field_t vtk_data[num_vtk_data]; /**< Data fields for the forest vtk writer pointing to the data above. */
};

/** Read the process local .vtu file written by the inbuilt ASCII vtk writer.
 * \param [in] fileprefix The prefix of the files that was passed to the vtk writer.
 * \param [in] mpirank    The rank of this process.
 * \return The content of the file.
 */
static std::string
read_vtu_file (const std::string &fileprefix, const int mpirank)
{
  // The ASCII writer pads the rank in the file name with zeros to four digits.
  char filename[BUFSIZ];
  snprintf (filename, BUFSIZ, "%s_%04d.vtu", fileprefix.c_str (), mpirank);
  std::ifstream file (filename);
  EXPECT_TRUE (file.is_open ()) << "Could not open file " << filename;
  if (!file.is_open ()) {
    return "";
  }
  std::stringstream content;
  content << file.rdbuf ();
  return content.str ();
}

/** Convert \a field with \ref detail::create_vtk_buffer and detail::fill_vtk_data_field and check the result.
 * The converted field must contain exactly one value per written element. Values of \a field for elements that are
 * not written are dropped and written elements without a value in \a field are set to zero.
 * \tparam TMeshClass       The mesh handle class.
 * \tparam TValueType       The type of the values in the data field.
 * \param [in] mesh         The mesh the data field belongs to.
 * \param [in] field        The data field to convert.
 * \param [in] write_ghosts If true, the ghost elements are written, too.
 */
template <t8_mesh_handle::T8MeshType TMeshClass, t8_mesh_handle::T8VtkDataValue TValueType>
void
check_vtk_data_conversion (const TMeshClass &mesh, const t8_mesh_handle::vtk_data_field<TValueType> &field,
                           const bool write_ghosts)
{
  constexpr bool is_scalar = std::is_same_v<TValueType, double>;
  constexpr int num_components = is_scalar ? 1 : 3;
  // Test the detail functions
  std::vector<double> buffer = t8_mesh_handle::detail::create_vtk_buffer (mesh, field, write_ghosts);
  t8_vtk_data_field_t vtk_field;
  t8_mesh_handle::detail::fill_vtk_data_field (field, buffer, vtk_field);

  EXPECT_EQ (vtk_field.type, is_scalar ? T8_VTK_SCALAR : T8_VTK_VECTOR);
  EXPECT_STREQ (vtk_field.description, field.name.c_str ());
  EXPECT_EQ (vtk_field.data, buffer.data ());
  const t8_locidx_t num_written_elements
    = write_ghosts ? mesh.get_num_local_elements () + mesh.get_num_ghosts () : mesh.get_num_local_elements ();
  ASSERT_EQ (buffer.size (), static_cast<size_t> (num_written_elements * num_components));
  for (t8_locidx_t ielem = 0; ielem < num_written_elements; ++ielem) {
    for (int icomponent = 0; icomponent < num_components; ++icomponent) {
      double expected_value = 0.0;
      if (static_cast<size_t> (ielem) < field.values.size ()) {
        if constexpr (is_scalar) {
          expected_value = field.values[ielem];
        }
        else {
          expected_value = field.values[ielem][icomponent];
        }
      }
      EXPECT_EQ (buffer[ielem * num_components + icomponent], expected_value) << "ielem = " << ielem;
    }
  }
}

/** Check that the mesh can be written with the standard output function \ref t8_mesh_handle::write_mesh_to_vtk. */
TEST_P (t8_mesh_io_test, write_mesh_to_vtk)
{
  EXPECT_TRUE (t8_mesh_handle::write_mesh_to_vtk (*mesh, output_prefix ("mesh_io")));
}

/** Check the conversion of the data fields of the handle into the data fields expected by the forest vtk writer.
 * The values of a data field can be given for the local elements only or for the local and the ghost elements.
 */
TEST_P (t8_mesh_io_test, to_vtk_data_field)
{
  const std::span<const double> local_volumes (volumes.data (), num_local_elements);
  const std::span<const t8_3D_vec> local_centroids (centroids.data (), num_local_elements);
  for (const bool write_ghosts : { false, true }) {
    // --- Scalar data fields. ---
    check_vtk_data_conversion (*mesh, t8_mesh_handle::vtk_scalar_field { "volume", volumes }, write_ghosts);
    check_vtk_data_conversion (*mesh, t8_mesh_handle::vtk_scalar_field { "local volume", local_volumes }, write_ghosts);
    // --- Vector data fields. ---
    check_vtk_data_conversion (*mesh, t8_mesh_handle::vtk_vector_field { "centroid", centroids }, write_ghosts);
    check_vtk_data_conversion (*mesh, t8_mesh_handle::vtk_vector_field { "local centroid", local_centroids },
                               write_ghosts);
  }
  // The longest name that fits into the description of the forest data field must be accepted.
  // Longer names abort with SC_CHECK_ABORTF, which cannot be tested here since it aborts all MPI processes.
  const std::string longest_name (T8_VTK_MAX_STRING_LENGTH - 1, 'a');
  check_vtk_data_conversion (*mesh, t8_mesh_handle::vtk_scalar_field { longest_name, volumes }, false);
}

/** Check the extended output function \ref t8_mesh_handle::write_mesh_to_vtk for different options and ghost element handling. */
TEST_P (t8_mesh_io_test, write_mesh_to_vtk_with_options)
{
  const std::span<const double> local_volumes (volumes.data (), num_local_elements);
  // Default options and no data fields.
  EXPECT_TRUE (t8_mesh_handle::write_mesh_to_vtk (*mesh, output_prefix ("with_options"), {}, {}, {}));

  for (const bool write_ghosts : { false, true }) {
    for (const bool do_not_use_API : { false, true }) {
      // All other options keep their default values.
      const t8_mesh_handle::vtk_write_options options { .write_ghosts = write_ghosts,
                                                        .do_not_use_API = do_not_use_API };

      // Only a scalar data field with values for the local elements only.
      EXPECT_TRUE (t8_mesh_handle::write_mesh_to_vtk (*mesh, output_prefix ("with_options"),
                                                      { { "volume", local_volumes } }, {}, options))
        << "write_ghosts = " << write_ghosts << ", do_not_use_API = " << do_not_use_API;
      // Scalar and vector data fields with values for the local and the ghost elements.
      EXPECT_TRUE (t8_mesh_handle::write_mesh_to_vtk (*mesh, output_prefix ("with_options"),
                                                      { { "level", levels }, { "volume", volumes } },
                                                      { { "centroid", centroids } }, options))
        << "write_ghosts = " << write_ghosts << ", do_not_use_API = " << do_not_use_API;
    }
  }
}

/** Check that the output functions of the handle write the same files as the forest output function
 * \ref t8_forest_write_vtk_ext with manually filled data fields.
 * Therefore, the process local files of the inbuilt ASCII writer are compared, since its output does not depend on
 * whether VTK is linked or not.
 */
TEST_P (t8_mesh_io_test, compare_output_to_forest)
{
  const std::string prefix_forest = "t8_gtest_mesh_io_forest";
  const std::string prefix_ext_data = "t8_gtest_mesh_io_handle_ext_data";
  // Use different values for the options to detect if they are mixed up.
  const bool write_treeid = true;
  const bool write_mpirank = false;
  const bool write_level = true;
  const bool write_element_id = false;
  const bool write_curved = false;
  const bool do_not_use_API = true;
  for (const bool write_ghosts : { false, true }) {
    // --- Write with the forest directly. ---
    EXPECT_TRUE (t8_forest_write_vtk_ext (mesh->get_forest (), prefix_forest.c_str (), write_treeid, write_mpirank,
                                          write_level, write_element_id, write_ghosts, write_curved, do_not_use_API,
                                          num_vtk_data, vtk_data));

    const t8_mesh_handle::vtk_write_options options { .write_ghosts = write_ghosts };
    // --- Write with the mesh handle. ---
    EXPECT_TRUE (t8_mesh_handle::write_mesh_to_vtk (*mesh, prefix_ext_data.c_str (),
                                                    { { "level", levels }, { "volume", volumes } },
                                                    { { "centroid", centroids } }, options));
    // --- Compare the process local files. ---
    const std::string forest_output = read_vtu_file (prefix_forest, mpirank);
    EXPECT_FALSE (forest_output.empty ());
    EXPECT_EQ (forest_output, read_vtu_file (prefix_ext_data, mpirank)) << "write_ghosts = " << write_ghosts;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_mesh_io, t8_mesh_io_test, testing::Combine (AllEclasses, testing::Range (2, 3)));

/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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
 * \file t8_gtest_custom_competence.cxx
 * Checks that custom-defined competences can be used as template parameters for the mesh. 
 */

#include <gtest/gtest.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8.h>

#include <t8_mesh_handle/mesh.hxx>
#include <t8_mesh_handle/element.hxx>
#include <t8_mesh_handle/competences.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_types/t8_operators.hxx>

/** Custom competence that needs access to the public members of the elements. */
template <typename TUnderlying>
struct dummy_get_level: public t8_crtp_operator<TUnderlying, dummy_get_level>
{
 public:
  t8_element_level
  get_level_dummy () const
  {
    auto forest = this->underlying ().get_mesh ()->get_forest ();
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, this->underlying ().get_tree_id ());
    const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, this->underlying ().get_tree_id (),
                                                                      this->underlying ().get_element_id ());
    return t8_forest_get_scheme (forest)->element_get_level (tree_class, element);
  }
};

/** Second custom competence. */
template <typename TUnderlying>
struct dummy_trivial: public t8_crtp_operator<TUnderlying, dummy_trivial>
{
 public:
  int
  get_value_dummy () const
  {
    return 1;
  }
};

/** This tests checks that custom defined competences can be used for the mesh class 
 *  and that we can use the functionality defined in the competence. 
 * Also checks that we can use more than one custom competence and that predefined competences can be additionally used.
 */
TEST (t8_gtest_custom_competence, custom_competence)
{
  // Define forest to construct mesh.
  const int level = 1;
  t8_cmesh_t cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  ASSERT_EQ (true, t8_forest_is_committed (forest));

  // Check mesh with custom defined competence.
  t8_mesh_handle::mesh<dummy_get_level> mesh = t8_mesh_handle::mesh<dummy_get_level> (forest);

  for (auto it = mesh.begin (); it != mesh.end (); ++it) {
    EXPECT_EQ (it->get_level (), it->get_level_dummy ());
    EXPECT_EQ (level, it->get_level_dummy ());
  }

  // Test with two custom competences and a predefined competence.
  using mesh_class = t8_mesh_handle::mesh<dummy_get_level, dummy_trivial, t8_mesh_handle::cache_centroid>;
  mesh_class mesh_more_competences = mesh_class (forest);

  for (auto it = mesh_more_competences.begin (); it != mesh_more_competences.end (); ++it) {
    EXPECT_EQ (it->get_level (), it->get_level_dummy ());
    EXPECT_EQ (it->get_value_dummy (), 1);
    EXPECT_FALSE (it->centroid_cache_filled ());
    auto centroid = it->get_centroid ();
    for (const auto &coordinate : centroid) {
      EXPECT_GE (1, coordinate);
      EXPECT_LE (0, coordinate);
    }
    EXPECT_TRUE (it->centroid_cache_filled ());
  }

  // Unref the forest.
  t8_forest_unref (&forest);
}

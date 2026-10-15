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
 * \file TODO
 *
 */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_forest/t8_forest_mesh_quality/t8_forest_aspect_ratio.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_schemes/t8_default/t8_default.hxx>

/** TODO
 */
TEST (t8_gtest_mesh_quality, aspect_ratio)
{
  // I wrote a dummy test only to verify that you can cll this function in general here. Of course the test does not make sense.
  const int level = 2;
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);
  t8_cmesh_new_from_class (cmesh, T8_ECLASS_QUAD, sc_MPI_COMM_WORLD);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);

  EXPECT_NEAR (0., t8_forest_average_aspect_ratio (forest), T8_PRECISION_SQRT_EPS);

  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); ++itree) {
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ++ielem) {
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      const auto elem_volume = t8_forest_element_volume (forest, itree, elem);
      // This of course only works because the aspect ratio currently returns the volume as a dummy.
      EXPECT_NEAR (elem_volume, t8_forest_element_aspect_ratio (forest, itree, elem), T8_PRECISION_SQRT_EPS);
    }
  }
  t8_forest_unref (&forest);
}

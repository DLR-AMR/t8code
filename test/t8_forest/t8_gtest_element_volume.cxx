/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

#include <sc_functions.h>
#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <test/t8_gtest_schemes.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_schemes/t8_default/t8_default_pyramid/t8_dpyramid_bits.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <test/t8_gtest_macros.hxx>

/**
 * This file tests the volume-computation of elements.
 */

/* Construct a forest of a hypercube with volume 1. If the element are refined uniformly
 * all elements have volume 1/global_num_elements. */

class t8_forest_volume: public testing::TestWithParam<std::tuple<std::tuple<int, t8_eclass_t>, int>> {
 protected:
  void
  SetUp () override
  {
    const int scheme_id = std::get<0> (std::get<0> (GetParam ()));
    scheme = create_from_scheme_id (scheme_id);
    eclass = std::get<1> (std::get<0> (GetParam ()));
    level = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (eclass, sc_MPI_COMM_WORLD, 0, 0, 0);
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, sc_MPI_COMM_WORLD);
  }
  void
  TearDown () override
  {
    t8_forest_unref (&forest);
  }
  t8_forest_t forest;
  const t8_scheme *scheme;
  t8_eclass_t eclass;
  int level;
};

/**
 * Compute the volume of a pyramid descending of a root-pyramid with volume 1/3
 * Pyramids need a special handling of the control-volume computation, because
 * they subdivide into pyramids and tetrahedra. Therefore in every refinement three 
 * types of elements occur:
 * 
 * 1. A pyramid with 1/8 of its parents volume
 * 2. A tetrahedron with a pyramid parent, having 1/16th of its parents volume.
 * 3. A tetrahedron with a tet-parent, having 1/8th of its parents volume.
 * 
 * On a leaf-level we therefore can have many different volumes for the
 * elements and compute it element-specific. 
 * \param[in] pyra A pyramid
 * \return The volume of the pyramid 
 */
double
pyramid_control_volume (t8_dpyramid_t *pyra)
{
  double control_volume = 1.0 / 3.0;
  /* Both pyramids and tets have 1/8th of the parents volume, if the shape does not switch. */
  control_volume /= 1 << ((pyra->pyramid.level) * 3);
  /* Ancestors switch the shape. A tetrahedron has a 1/16th of its parents volume. 
   * For all levels we already divided the control-volume by 8, hence we 
   * divide it by 2 once. */
  if (pyra->switch_shape_at_level > 0) {
    control_volume /= 2;
  }

  return control_volume;
}

TEST_P (t8_forest_volume, volume_check)
{
  /* Compute the global number of elements */
  const t8_gloidx_t global_num_elements = t8_forest_get_global_num_elements (forest);
  /* Vertices have a volume of 0. */
  const double control_volume = (eclass == T8_ECLASS_VERTEX) ? 0.0 : (1.0 / global_num_elements);

  ASSERT_EQ (t8_forest_get_dimension (forest), t8_cmesh_get_dimension (t8_forest_get_cmesh (forest)));

  const t8_locidx_t local_num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over all elements. */
  for (t8_locidx_t itree = 0; itree < local_num_trees; itree++) {
    const t8_locidx_t tree_elements = t8_forest_get_tree_num_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < tree_elements; ielement++) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      const double volume = t8_forest_element_volume (forest, itree, element);
      if (eclass == T8_ECLASS_PYRAMID) {
        const double shape_volume = pyramid_control_volume ((t8_dpyramid_t *) element);
        EXPECT_NEAR (volume, shape_volume, T8_PRECISION_SQRT_EPS);
      }
      else {
        EXPECT_NEAR (volume, control_volume, T8_PRECISION_SQRT_EPS);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_element_volume, t8_forest_volume,
                          testing::Combine (AllSchemes, testing::Range (0, 4)));

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

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>

/* In this file we collect tests for t8code's cmesh geometry module.
 * These tests are
 *  - geometry.geometry_linear:  Check that linear geometry has correct name and dimension.
 *  - geometry.geometry_zero:    Check that zero geometry has correct name and dimension. 
 *  - geometry.cmesh_geometry_linear: For each dimension create a cmesh with linear geometry.
 *                                   We then create random points and check whether their geometry
 *                                   is computed correctly.
 *  - test_geometry.cmesh_geometry: 
 *  - test_geometry.cmesh_geometry_unique: Check that we can acces the geometry via the tree id if
 *                                   we only use one geometry and did not specify tree ids for it.
 *                                   In this case t8code should automatically associate this geometry to all trees.
 *  - test_geometry.geom_handler_register: Tests the geometry_handler register and find interface.
 */
 /* TODO: 
  * - Add a test for the jacobian, as soon as its implemented in parameterized test geometry.cmesh_geometry_linear.
  * - We need to use a different geometry since we remove the identity geomeytry in test_cmesh_geometry. 
  * - We need to use a different geometry since we remove the identity geomeytry in test_cmesh_geometry_unique.
  */

/* Check that the linear geometry for dimensions 0,1,2,3
 * has the correct name and dimension. */

/* *INDENT-OFF* */
class geometry:public testing::TestWithParam <int> {
protected:
  void SetUp () override {
    dim = GetParam();
  }
  int   dim;

};
/* *INDENT-ON* */

/* Check that the linear geometry for dimensions 0,1,2,3
 * has the correct name and dimension. */
TEST_P (geometry, geometry_linear)
{
  t8_geometry_linear  linear_geom (dim);
  char                name[BUFSIZ];
  snprintf (name, BUFSIZ, "t8_geom_linear_%i", dim);
  EXPECT_FALSE (strcmp (linear_geom.t8_geom_get_name (), name));
  EXPECT_EQ (dim, linear_geom.t8_geom_get_dimension ());
}

/* Check that the zero geometry for dimensions 0,1,2,3
 * has the correct name and dimension. */
TEST_P (geometry, geometry_zero)
{
  t8_geometry_zero    zero_geom (dim);
  char                name[BUFSIZ];
  snprintf (name, BUFSIZ, "t8_geom_zero_%i", dim);
  EXPECT_FALSE (strcmp (zero_geom.t8_geom_get_name (), name));
  EXPECT_EQ (dim, zero_geom.t8_geom_get_dimension ());
}

/* Check whether the linear geometry map is correct.
 * We create a cmesh of one tree and unit point/line/square/cube
 * geometry. We then create random points in a reference tree and
 * check whether the evaluation is correct. */
TEST_P (geometry, cmesh_geometry_linear)
{
  /* TODO: Add a test for the jacobian, as soon as its implemented. */
  /* Create random points in [0,1]^d and check if they are mapped correctly. */
  t8_geometry_linear  linear_geom (dim);
  t8_cmesh_t          cmesh;

  int                 num_points = 10000;
  double              point[3];
  double              point_mapped[3];
  int                 seed = 0; /* RNG seed */
  t8_eclass_t         tree_class;
  const t8_geometry_c *cmesh_geom;

  /* Build a one tree cmesh on the unit square with linear geometry. */
  switch (dim) {
  case 0:
    tree_class = T8_ECLASS_VERTEX;
    break;
  case 1:
    tree_class = T8_ECLASS_LINE;
    break;
  case 2:
    tree_class = T8_ECLASS_QUAD;
    break;
  case 3:
    tree_class = T8_ECLASS_HEX;
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
  cmesh = t8_cmesh_new_hypercube (tree_class, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Double check that the geometry is the linear geometry. */
  cmesh_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  EXPECT_TRUE (cmesh_geom != NULL);
  EXPECT_FALSE (strcmp (cmesh_geom->t8_geom_get_name (),
                        linear_geom.t8_geom_get_name ()));

  srand (seed);
  for (int ipoint = 0; ipoint < num_points; ++ipoint) {
    /* Compute random coordinates in [0,1].
     * These are seen as reference coordinates in the single
     * cmesh tree. Our geometry will map them into the physical
     * space. Since this space is also [0,1] and the cmesh only
     * has one tree, the mapped coordinates must be the same as the 
     * reference coordinates. */
    point[0] = (double) rand () / RAND_MAX;
    point[1] = (double) rand () / RAND_MAX;
    point[2] = (double) rand () / RAND_MAX;

    /* Evaluate the geometry */
    t8_geometry_evaluate (cmesh, 0, point, point_mapped);
    /* Check that the first dim coordinates are the same */
    int                 idim;
    for (idim = 0; idim < dim; ++idim) {
      const double        tolerance = 1e-14;
      EXPECT_TRUE (fabs (point[idim] - point_mapped[idim]) < tolerance);
    }
    /* Check that the remaining entries are 0. */
    for (; idim < 3; ++idim) {
      EXPECT_EQ (point_mapped[idim], 0);
    }
  }
  /* Destroy the cmesh */
  t8_cmesh_destroy (&cmesh);
}

/* *INDENT-OFF* */
INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry, geometry,testing::Range(0,4));
/* *INDENT-ON* */

TEST (test_geometry, cmesh_geometry)
{
  t8_cmesh_t          cmesh;
  /* TODO: We need to use a different geometry since we remove the identity geomeytry. */

  t8_geometry_linear *linear_geom = new t8_geometry_linear (2);
  t8_geometry_zero   *zero_geom = new t8_geometry_zero (2);
  const t8_geometry_c *found_geom;

  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  /* Register the linear geometry and zero geometry to this cmesh. */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  t8_cmesh_register_geometry (cmesh, zero_geom);
  /* Set the id geometry for the trees. */
  t8_cmesh_set_tree_geometry (cmesh, 0, linear_geom->t8_geom_get_name ());
  t8_cmesh_set_tree_geometry (cmesh, 1, zero_geom->t8_geom_get_name ());
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id. */
  found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  EXPECT_TRUE (found_geom != NULL);
  EXPECT_TRUE (found_geom == linear_geom);

  found_geom = t8_cmesh_get_tree_geometry (cmesh, 1);
  EXPECT_TRUE (found_geom != NULL);
  EXPECT_TRUE (found_geom == zero_geom);
  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, cmesh_geometry_unique)
{
  t8_cmesh_t          cmesh;
  /* TODO: We need to use a different geometry since we remove the identity geomeytry. */

  t8_geometry_linear *linear_geom = new t8_geometry_linear (2);
  const t8_geometry_c *found_geom;

  /* Build a simple 2 tree cmesh and set geometries for the trees. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  /* Register the linear_geometry to this cmesh. */
  t8_cmesh_register_geometry (cmesh, linear_geom);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);

  /* Check that we can get the geometry back over the tree id.
   * This must now work even though we did not register the geometry for 
   * this tree. Since we only have one geometry. */
  found_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  EXPECT_TRUE (found_geom != NULL);
  EXPECT_EQ (found_geom, linear_geom);

  /* clean-up */
  t8_cmesh_destroy (&cmesh);
}

TEST (test_geometry, geom_handler_register)
{
  t8_geometry_handler_t *geom_handler;
  const t8_geometry_c *found_geom;

  /* Initialize a geometry handler. */
  t8_geom_handler_init (&geom_handler);

  /* For each dimension build the zero geometry and register it.
   * We then commit the handler and check that we can find the geometries. */
  for (int idim = 0; idim <= 3; ++idim) {
    t8_geometry_zero   *zero_geom = new t8_geometry_zero (idim);
    /* Register the geometry. */
    t8_geom_handler_register_geometry (geom_handler, zero_geom);
  }
  /* Commit the handler */
  t8_geom_handler_commit (geom_handler);

  /* Check find geometry. */
  for (int idim = 0; idim < 3; ++idim) {
    t8_geometry_zero    zero_geom (idim);
    const char         *name;

    /* Get the name of this geometry. */
    name = zero_geom.t8_geom_get_name ();

    /* Find the geometry by name. */
    found_geom = t8_geom_handler_find_geometry (geom_handler, name);
    EXPECT_TRUE (found_geom != NULL);
    EXPECT_FALSE (strcmp (found_geom->t8_geom_get_name (), name));
  }
  /* Try to find a different geometry. Must return NULL. */
  found_geom =
    t8_geom_handler_find_geometry (geom_handler, "random_name34823412414");
  EXPECT_TRUE (found_geom == NULL);

  /* clean-up */
  t8_geom_handler_destroy (&geom_handler);
  EXPECT_TRUE (geom_handler == NULL);
}

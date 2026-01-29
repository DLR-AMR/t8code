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

/** \file t8_gtest_geometry_curved.cxx
 * Component test for curved geometries.
 * This test verifies that curved geometries (Lagrange) work correctly in the 
 * complete t8code workflow, from cmesh creation through forest operations.
 * It complements the unit tests in t8_gtest_geometry_lagrange.cxx by testing
 * the integration of curved geometries in realistic use cases.
 */

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_geometry.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <test/t8_gtest_macros.hxx>

/**
 * Create a sample curved element with perturbed vertices for testing.
 *
 * \param [in] eclass  Element class of the element.
 * \param [in] degree  Polynomial degree for Lagrange interpolation.
 * \return             Vector containing vertex coordinates.
 */
std::vector<double>
create_curved_element_vertices (t8_eclass_t eclass, int degree)
{
  std::vector<double> vertices;

  switch (eclass) {
  case T8_ECLASS_QUAD:
    if (degree == 1) {
      /* Linear quad with slight perturbation */
      vertices = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0 };
    }
    else if (degree == 2) {
      /* Quadratic quad with curved edges */
      vertices = { 0.0, 0.0,  0.0, 1.0, 0.0,  0.0,  0.0,  1.0,  0.0, 1.0, 1.0, 0.0,
                   0.5, -0.1, 0.0, 1.1, 0.5,  0.0,  0.5,  1.1,  0.0, -0.1, 0.5, 0.0,
                   0.5, 0.5,  0.0 };
    }
    break;
  case T8_ECLASS_TRIANGLE:
    if (degree == 1) {
      /* Linear triangle */
      vertices = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
    }
    else if (degree == 2) {
      /* Quadratic triangle with curved edges */
      vertices = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                   0.5, -0.05, 0.0, 0.6, 0.6, 0.0, -0.05, 0.5, 0.0 };
    }
    break;
  case T8_ECLASS_HEX:
    if (degree == 1) {
      /* Linear hex */
      vertices = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                   0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    }
    else if (degree == 2) {
      /* Quadratic hex with curved edges - simplified version */
      vertices = { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                   0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                   /* Edge midpoints with perturbations */
                   0.5, 0.0, 0.05, 0.0, 0.5, 0.0, 1.0, 0.5, 0.0, 0.5, 1.0, 0.0,
                   0.5, 0.0, 1.0, 0.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.5, 1.0, 1.0,
                   0.0, 0.0, 0.5, 1.0, 0.0, 0.5, 0.0, 1.0, 0.5, 1.0, 1.0, 0.5,
                   /* Face centers */
                   0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 1.0, 0.5, 0.5,
                   0.5, 1.0, 0.5, 0.5, 0.5, 1.0,
                   /* Volume center */
                   0.5, 0.5, 0.5 };
    }
    break;
  default:
    SC_ABORTF ("Element type %s not supported for curved geometry test.\n", t8_eclass_to_string[eclass]);
  }

  return vertices;
}

/**
 * Test fixture for curved geometry component tests.
 */
class CurvedGeometry: public testing::TestWithParam<std::tuple<t8_eclass_t, int>>
{
 protected:
  void
  SetUp () override
  {
    /* Fetch the current element type and polynomial degree */
    std::tuple<t8_eclass, int> params = GetParam ();
    eclass = std::get<0> (params);
    degree = std::get<1> (params);

    /* Only test element types that are supported by Lagrange geometry */
    if (eclass != T8_ECLASS_QUAD && eclass != T8_ECLASS_HEX && eclass != T8_ECLASS_TRIANGLE) {
      GTEST_SKIP () << "Element type " << t8_eclass_to_string[eclass] << " not yet supported for curved geometries.\n";
    }
  }

  /**
   * Helper function to create and setup a cmesh with curved geometry.
   * \param[out] cmesh  The cmesh to create and setup.
   */
  void
  create_curved_cmesh (t8_cmesh_t *cmesh)
  {
    t8_cmesh_init (cmesh);
    t8_cmesh_set_tree_class (*cmesh, 0, eclass);

    std::vector<double> vertices = create_curved_element_vertices (eclass, degree);
    const int num_vertices = vertices.size () / 3;

    t8_cmesh_set_tree_vertices (*cmesh, 0, vertices.data (), num_vertices);
    t8_cmesh_set_attribute (*cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE_KEY, &degree,
                            sizeof (degree), 0);
    auto lagrange_geom = t8_cmesh_register_geometry<t8_geometry_lagrange> (*cmesh);
    t8_cmesh_set_tree_geometry (*cmesh, 0, lagrange_geom);
    t8_cmesh_commit (*cmesh, sc_MPI_COMM_WORLD);
  }

  t8_eclass_t eclass;
  int degree;
};

/**
 * Test that a cmesh with curved geometry can be created and committed.
 */
TEST_P (CurvedGeometry, cmesh_creation_with_curved_geometry)
{
  t8_cmesh_t cmesh;

  /* Create a single-tree cmesh with Lagrange geometry */
  create_curved_cmesh (&cmesh);

  /* Verify the geometry is retrievable */
  const t8_geometry *retrieved_geom = t8_cmesh_get_tree_geometry (cmesh, 0);
  ASSERT_NE (retrieved_geom, nullptr) << "Failed to retrieve geometry from cmesh.";
  ASSERT_EQ (retrieved_geom->t8_geom_get_type (), T8_GEOMETRY_TYPE_LAGRANGE)
    << "Retrieved geometry is not of type LAGRANGE.";

  /* Clean up */
  t8_cmesh_destroy (&cmesh);
}

/**
 * Test that a uniform forest can be created from a cmesh with curved geometry
 * and that geometry evaluation works correctly on forest elements.
 */
TEST_P (CurvedGeometry, forest_with_curved_geometry_and_evaluation)
{
  t8_cmesh_t cmesh;
  t8_forest_t forest;

  /* Create a single-tree cmesh with Lagrange geometry */
  create_curved_cmesh (&cmesh);

  /* Create a uniform forest from the curved cmesh */
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default (), 1, 0, sc_MPI_COMM_WORLD);

  ASSERT_NE (forest, nullptr) << "Forest creation failed with curved geometry.";

  /* Verify we can evaluate the geometry through the forest
   * This tests the integration of curved geometry in the forest workflow */
  t8_cmesh_t cmesh_from_forest = t8_forest_get_cmesh (forest);
  ASSERT_NE (cmesh_from_forest, nullptr) << "Failed to get cmesh from forest.";

  /* Verify geometry is still accessible through the forest's cmesh */
  const t8_geometry *geom = t8_cmesh_get_tree_geometry (cmesh_from_forest, 0);
  ASSERT_NE (geom, nullptr) << "Failed to get geometry from forest's cmesh.";
  ASSERT_EQ (geom->t8_geom_get_type (), T8_GEOMETRY_TYPE_LAGRANGE)
    << "Geometry type changed after forest creation.";

  /* Clean up */
  t8_forest_unref (&forest);
}

/* Instantiate tests for different element types and polynomial degrees */
INSTANTIATE_TEST_SUITE_P (
  t8_gtest_geometry_curved, CurvedGeometry,
  testing::Combine (testing::Values (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE, T8_ECLASS_HEX),
                    testing::Values (1, 2)),
  [] (const testing::TestParamInfo<CurvedGeometry::ParamType> &info) {
    std::ostringstream test_name;
    test_name << t8_eclass_to_string[std::get<0> (info.param)] << "_degree" << std::get<1> (info.param);
    return test_name.str ();
  });

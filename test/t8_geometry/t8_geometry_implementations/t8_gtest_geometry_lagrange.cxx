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

/** \file t8_gtest_geometry_lagrange.cxx
* Provide tests to check the geometry mappings via Lagrange interpolation.
* The tests are parameterized on the element type and the polynomial degree.
*/

#include <array>
#include <cstdlib>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>
#include <test/t8_gtest_macros.hxx>
#include <t8_eclass.h>
#include <t8_types/t8_vec.hxx>
#include <t8_cmesh.hxx>
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>

/**
 * Generate a random double precision number in [0, 1].
 *
 * \return  Random number.
 */
double
random_number ()
{
  return static_cast<double> (std::rand ()) / RAND_MAX;
}

/**
 * Add random perturbation to a vector.
 *
 * \param vec            Vector to be perturbed.
 * \param max_amplitude  Value of the maximum perturbation.
 * \return               Perturbed vector.
 */
std::vector<double>
perturb (const std::vector<double> &vec, const double max_amplitude)
{
  std::vector<double> perturbed_vec;
  perturbed_vec.reserve (vec.size ());
  for (const auto &v : vec)
    perturbed_vec.push_back (v + random_number () * max_amplitude);
  return perturbed_vec;
}

/**
 * Return true if two arrays are element-wise equal within a tolerance.
 * \tparam U   Iterable container with size.
 * \tparam V   Iterable container with size.
 * \param a    First array.
 * \param b    Second arary.
 * \param tol  Absolute tolerance.
 * \return     True if all the elements of \a a and \a b are close to each
 *             other, false otherwise.
 */
template <typename U, typename V>
bool
allclose (const U &a, const V &b, double tol = T8_PRECISION_SQRT_EPS)
{
  T8_ASSERT (a.size () == b.size ());
  for (size_t i = 0; i < a.size (); ++i)
    if (fabs (a[i] - b[i]) > tol)
      return false;
  return true;
}

/**
 * Sample random points in the reference domain of an element class.
 *
 * \param eclass   Element class.
 * \param n_point  Number of points to generate.
 * \return         Coordinates of the points, given in x,y,z.
 */
std::vector<std::array<double, T8_ECLASS_MAX_DIM>>
sample (t8_eclass_t eclass, uint32_t n_point)
{
  std::srand (time (NULL));
  std::vector<std::array<double, T8_ECLASS_MAX_DIM>> points (n_point);
  switch (eclass) {
  case T8_ECLASS_LINE:
    for (auto &pt : points) {
      pt[0] = random_number ();
      pt[1] = 0;
      pt[2] = 0;
    }
    break;
  case T8_ECLASS_QUAD:
    for (auto &pt : points) {
      pt[0] = random_number ();
      pt[1] = random_number ();
      pt[2] = 0;
    }
    break;
  case T8_ECLASS_HEX:
    for (auto &pt : points) {
      pt[0] = random_number ();
      pt[1] = random_number ();
      pt[2] = random_number ();
    }
    break;
  default:
    SC_ABORTF ("Sampling on element class %d not supported.\n", eclass);
  }
  return points;
}

/**
 * Create a sample t8_lagrange_element.
 *
 * The goal of this function is to quickly instantiate a t8_lagrange_element
 * for the purpose of testing.
 *
 * \param eclass  Element class of the element.
 * \param degree  Polynomial degree.
 * \return        t8_lagrange_element.
 */
t8_lagrange_element
create_sample_element (t8_eclass_t eclass, int degree)
{
  std::ostringstream invalid_degree;
  invalid_degree << "Degree " << degree << " is not yet supported for " << t8_eclass_to_string[eclass]
                 << " elements.\n";
  std::vector<double> vertices;
  switch (eclass) {
  case T8_ECLASS_TRIANGLE:
    switch (degree) {
    case 1:
      vertices = { 0, 0, 0, 1, 0, 0, 1, 1, 0 };
      break;
    case 2:
      vertices = { 5.0, 0.0, 0.0, 2.0, 3.0, 0.0, 1.0, 1.0, 0.0, 3.5, 1.5, 0, 1.5, 2, 0, 3, -0.5, 0 };
      break;
    default:
      SC_ABORT (invalid_degree.str ().c_str ());
    }
    break;
  case T8_ECLASS_QUAD:
    switch (degree) {
    case 1:
      vertices = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.8, 1, 0 };
      break;
    case 2:
      vertices
        = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.8, 1, 0, 0, 0.5, 0, 0.95, 0.5, 0, 0.5, 0, 0, 0.4, 1.1, 0, 0.45, 0.5, 0 };
      break;
    default:
      SC_ABORT (invalid_degree.str ().c_str ());
    }
    break;
  case T8_ECLASS_HEX:
    switch (degree) {
    case 1:
      vertices = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
      vertices = perturb (vertices, 1 / 5.0 * (1 / degree));
      break;
    case 2:
      /* clang-format off */
      vertices = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1,   0, 0.5, 1, 1, 0.5, 1, 0.5, 0, 1, 0.5, 1, 1, 0.5, 0.5, 0.5, 0, 0, 0.5, 0, 1, 0.5, 0, 0.5, 0.5, 1, 0, 0.5, 1, 1, 0.5, 1, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 1, 0.5, 0.5, 0.5 };
      /* clang-format off */
      vertices = perturb (vertices, 1 / 5.0 * (1 / degree));
      break;
    default:
      SC_ABORT (invalid_degree.str ().c_str ());
    }
    break;
  default:
    SC_ABORTF ("Not implemented for %s elements.\n", t8_eclass_to_string[eclass]);
  }
  return t8_lagrange_element (eclass, degree, vertices);
}

/**
 * Common resources for all the tests.
 *
 */
class LagrangeCmesh: public testing::TestWithParam<std::tuple<t8_eclass_t, int>> {
 protected:
  void
  SetUp () override
  {
    /* Fetch the current element type and polynomial degree */
    std::tuple<t8_eclass, int> params = GetParam ();
    eclass = std::get<0> (params);
    degree = std::get<1> (params);
    if (eclass != T8_ECLASS_TRIANGLE && eclass != T8_ECLASS_QUAD && eclass != T8_ECLASS_HEX)
      GTEST_SKIP () << "Element type not yet implemented.\n";
  }

  t8_eclass_t eclass;
  int degree;
};

/**
 * Main test to check the correctness of the Lagrange geometries.
 *
 */
TEST_P (LagrangeCmesh, lagrange_mapping)
{
  /* Create a coarse mesh consisting of one single element
   * (one element is sufficient: we want to test the mapping) */
  t8_lagrange_element lag = create_sample_element (eclass, degree);
  /* Compare the mappings on each boundary face of the Lagrange element */
  std::vector<t8_lagrange_element> faces = lag.decompose ();
  uint32_t i_face = 0;
  for (const auto &face : faces) {
    auto points_on_face = sample (face.get_type (), T8_NUM_SAMPLE_POINTS);
    for (const auto &point : points_on_face) {
      auto mapped1 = face.evaluate (point);
      auto mapped2 = lag.evaluate (face.map_on_face (lag.get_type (), i_face, point));
      ASSERT_TRUE (allclose (mapped1, mapped2));
    }
    ++i_face;
  }
}

/* clang-format off */
INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry_lagrange, LagrangeCmesh,
  testing::Combine (AllEclasses, testing::Range (1, T8_GEOMETRY_MAX_POLYNOMIAL_DEGREE + 1)),
  [] (const testing::TestParamInfo<LagrangeCmesh::ParamType> &info) {
    std::ostringstream test_name;
    test_name << t8_eclass_to_string[std::get<0> (info.param)]
              << "_degree" << std::get<1> (info.param);
    return test_name.str ();});
/* clang-format on */

#if T8_ENABLE_DEBUG

/**
 * Tests the compatibility checking for the Lagrange geometry.
 * The geometry should throw assertions if the geometry is not compatible with an assigned tree.
 */
TEST (test_geometry_lagrange, incompatible_geometry)
{
  t8_cmesh_t cmesh;
  int degree = 1;

  t8_debugf ("Testing geometry compatibility checking for lagrange geometry.\n");

  /* Build a simple set geometries for the tree. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_QUAD], 4);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE_KEY, &degree, sizeof (degree),
                          0);

  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  /* Register the t8_geometry_lagrange geometry to this cmesh. */
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh);
  /* Should return true since the t8_geometry_lagrange geometry is compatible with quads. */
  ASSERT_TRUE (t8_cmesh_validate_geometry (cmesh, 0));
  t8_cmesh_destroy (&cmesh);

  /* Build a simple set geometries for the tree. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_HEX], 8);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE_KEY, &degree, sizeof (degree),
                          0);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_PRISM);
  t8_cmesh_set_tree_vertices (cmesh, 1, *t8_element_corner_ref_coords[T8_ECLASS_PRISM], 6);
  t8_cmesh_set_attribute (cmesh, 1, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE_KEY, &degree, sizeof (degree),
                          0);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  /* Register the t8_geometry_lagrange to this cmesh.
   * We register it after committing because it would throw an assertion and we do not have death tests.*/
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh);
  /* Check validity after committing to circumvent the assertion.
   * Should return false since the t8_geometry_lagrange geometry is not compatible with prisms. */
  ASSERT_FALSE (t8_cmesh_validate_geometry (cmesh, 0));
  t8_cmesh_destroy (&cmesh);

  degree = T8_GEOMETRY_MAX_POLYNOMIAL_DEGREE + 1;
  /* Build a simple set geometries for the tree. */
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_vertices (cmesh, 0, *t8_element_corner_ref_coords[T8_ECLASS_HEX], 8);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE_KEY, &degree, sizeof (degree),
                          0);
  /* Commit the cmesh */
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  /* Register the t8_geometry_lagrange to this cmesh.
   * We register it after committing because it would throw an assertion and we do not have death tests.*/
  t8_cmesh_register_geometry<t8_geometry_lagrange> (cmesh);
  /* Check validity after committing to circumvent the assertion.
   * Should return false since the maximum polynomial degree is exceeded. */
  ASSERT_FALSE (t8_cmesh_validate_geometry (cmesh, 0));
  t8_cmesh_destroy (&cmesh);
}
#endif /* T8_ENABLE_DEBUG */

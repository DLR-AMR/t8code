/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2023 the developers

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
* The tests are parametrized on the element type and the polynomial degree.
*/

#include <sstream>
#include <vector>

#include <gtest/gtest.h>
#include <t8_eclass.h>
#include <t8_vec.h>
#include <t8_element_cxx.hxx>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_lagrange.hxx>
#include <test/t8_gtest_macros.hxx>

#if T8_ENABLE_LESS_TESTS
#define MAX_POLYNOMIAL_DEGREE 1
#else
#define MAX_POLYNOMIAL_DEGREE 1
#endif

/** Constructs a cmesh from a single Lagrange line element.
 * \param cmesh   Uninitialized cmesh.
 * \param degree  Polynomial degree.
 * \return        A committed t8_cmesh structure with a Lagrange geometry.
 */
t8_cmesh_t
t8_cmesh_new_lagrange_line (t8_cmesh_t cmesh, int degree)
{
  std::vector<double> vertices;
  switch (degree) {
  case 1:
    vertices = { 0, 0, 0, 1, 0, 0 };
    break;
  case 2:
    vertices = { -1, 0, 0, 1, 0, 0, 0, 0, 0 };
    break;
  default:
    SC_ABORT ("Invalid degree\n");
    return nullptr;
  }
  t8_geometry_c *geometry = new t8_geometry_lagrange (1);
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices.data (), vertices.size ());
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  return cmesh;
};

/** Constructs a cmesh from a single Lagrange triangle element.
 * \param cmesh   Uninitialized cmesh.
 * \param degree  Polynomial degree.
 * \return        A committed t8_cmesh structure with a Lagrange geometry.
 */
t8_cmesh_t
t8_cmesh_new_lagrange_tri (t8_cmesh_t cmesh, int degree)
{
  std::vector<double> vertices;
  switch (degree) {
  case 1:
    vertices = { 0, 0, 0, 1, 0, 0, 1, 1, 0 };
    break;
  case 2:
    vertices = { 5.0, 0.0, 0.0, 2.0, 3.0, 0.0, 1.0, 1.0, 0.0, 3.5, 1.5, 0, 1.5, 2, 0, 3, -0.5, 0 };
    break;
  default:
    SC_ABORT ("Invalid degree\n");
    return nullptr;
  }
  t8_geometry_c *geometry = new t8_geometry_lagrange (2);
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices.data (), vertices.size ());
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  return cmesh;
};

/** Constructs a cmesh from a single Lagrange quadrilateral element.
 * \param cmesh   Uninitialized cmesh.
 * \param degree  Polynomial degree.
 * \return        A committed t8_cmesh structure with a Lagrange geometry.
 */
t8_cmesh_t
t8_cmesh_new_lagrange_quad (t8_cmesh_t cmesh, int degree)
{
  std::vector<double> vertices;
  switch (degree) {
  case 1:
    vertices = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.8, 1, 0 };
    break;
  default:
    SC_ABORT ("Invalid degree\n");
    return nullptr;
  }
  t8_geometry_c *geometry = new t8_geometry_lagrange (2);
  t8_cmesh_init (&cmesh);
  t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
  t8_cmesh_register_geometry (cmesh, geometry);
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices.data (), vertices.size ());
  t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  return cmesh;
};

/** Constructs a cmesh from a single Lagrange element.
 * \param cmesh   Uninitialized cmesh.
 * \param eclass  Element type.
 * \param degree  Polynomial degree.
 * \return        A committed t8_cmesh structure with a Lagrange geometry.
 */
t8_cmesh_t
t8_cmesh_new_lagrange_elem (t8_cmesh_t cmesh, t8_eclass_t eclass, int degree)
{
  switch (eclass) {
  case T8_ECLASS_LINE:
    return t8_cmesh_new_lagrange_line (cmesh, degree);
    break;
  case T8_ECLASS_TRIANGLE:
    return t8_cmesh_new_lagrange_tri (cmesh, degree);
    break;
  case T8_ECLASS_QUAD:
    return t8_cmesh_new_lagrange_quad (cmesh, degree);
    break;
  default:
    SC_ABORTF ("Invalid element class %d.\n", eclass);
    return nullptr;
  }
};

/**
 * Probes points in the element and gives its expected mapped coordinates.
 * 
 * \param eclass  Element type.
 * \param degree  Polynomial degree.
 * \return        Sampling points and the mapped sampling points.
 */
std::tuple<std::vector<double>, std::vector<double>>
t8_cmesh_create_points (t8_eclass_t eclass, int degree)
{
  std::vector<double> ref_points;
  std::vector<double> mapped_points;
  std::ostringstream invalid_degree;
  invalid_degree << "Invalid degree " << degree << ".\n";
  switch (eclass) {
  case T8_ECLASS_LINE:
    if (degree == 1) {
      ref_points = { 5.5, 0, 0 };
      mapped_points = { 5.5, 0, 0 };
    }
    else
      SC_ABORTF (invalid_degree.str ().c_str ());
    break;
  case T8_ECLASS_TRIANGLE:
    if (degree == 1) {
      ref_points = { 0.5, 0.5, 0 };
      mapped_points = { 0.5, 0.5, 0 };
    }
    else if (degree == 2) {
      ref_points = { 0, 0, 0 };
      mapped_points = { 5.0, 0.0, 0.0, 2.0, 3.0, 0.0, 1.0, 1.0, 0.0, 3.5, 1.5, 0, 1.5, 2, 0, 3, -0.5, 0 };
    }
    else
      SC_ABORTF (invalid_degree.str ().c_str ());
    break;
  case T8_ECLASS_QUAD:
    if (degree == 1) {
      ref_points = { 1, 0.5, 0 };
      mapped_points = { 0.9, 0.5, 0 };
    }
    else
      SC_ABORTF (invalid_degree.str ().c_str ());
    break;
  default:
    SC_ABORTF ("Invalid element class %d.\n", eclass);
  }
  return std::make_tuple (ref_points, mapped_points);
};

class Boundary1D {
 public:
  /**
  * Construct a new Boundary1D object.
  * 
  * \param degree  Polynomial degree.
  */
  Boundary1D (uint degree): degree (degree)
  {
    /* Generate the equidistant nodes adhering to the numbering convention */
    switch (degree) {
    case 1:
      parametric_nodes = { 0, 0, 0, 1, 0, 0 };
      break;
    case 2:
      parametric_nodes = { 0, 0, 0, 1, 0, 0, 0.5, 0, 0 };
      break;
    default:
      SC_ABORTF ("Degree must be in [1, 2].\n");
    }
  };

  /**
   * Set the nodes of the boundary curve.
   * 
   * \param nodes  Coordinates in the real space, {x1, y1, z1, x2, ...}.
   */
  void
  setNodes (std::vector<double> &nodes)
  {
    if (nodes.size () != parametric_nodes.size ())
      SC_ABORTF ("Provide the 3 coordinates of the nodes.\n");
    /* Create a cmesh with a single line element */
    t8_geometry_c *geometry = new t8_geometry_lagrange (1);
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
    t8_cmesh_register_geometry (cmesh, geometry);
    t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_LINE);
    t8_cmesh_set_tree_vertices (cmesh, 0, nodes.data (), nodes.size ());
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  };

  /**
   * ...
   * 
   * \param  sampling_points  Parametric coordinates of the points to be mapped.
   * \return                  Mapped sampling points.
   */
  std::vector<double>
  sample (const std::vector<double> &sampling_points)
  {
    const size_t n_point = sampling_points.size ();
    std::vector<double> mapped (n_point * 3, 0);
    for (size_t i = 0; i < n_point; ++i) {
      double point = sampling_points[i];
      std::array<double, 3> ref_point = { point, 0, 0 };
      t8_geometry_evaluate (cmesh, 0, ref_point.data(), 1, mapped.data() + 3 * i);
    }
    return mapped;
  };

  ~Boundary1D ()
  {
    std::cout << "Destructor called.\n";
    t8_cmesh_destroy (&cmesh);
  }

 private:
  uint degree;
  std::vector<double> parametric_nodes;
  t8_cmesh_t cmesh;
};

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
    /* Create a coarse mesh consisting of one single element
     * (one element is sufficient: we want to test the mapping) */
    cmesh = t8_cmesh_new_lagrange_elem (cmesh, eclass, degree);
  }

  void
  TearDown () override
  {
    t8_cmesh_destroy (&cmesh);
  }

  t8_eclass_t eclass;
  int degree;
  t8_cmesh_t cmesh;
};

/**
 * Main test to check the correctness of the Lagrange geometries.
 * 
 */
TEST_P (LagrangeCmesh, lagrange_mapping)
{
  t8_errorf ("\n-------------------\nEclass: %s, degree: %d\n-------------------\n", t8_eclass_to_string[eclass],
             degree);
  std::tuple<std::vector<double>, std::vector<double>> res = t8_cmesh_create_points (eclass, degree);
  const auto ref_points = std::get<0> (res);
  const auto expected = std::get<1> (res);
  double mapped_points[3];
  t8_geometry_evaluate (cmesh, 0, ref_points.data (), 1, mapped_points);
  for (size_t coord = 0; coord < 3; ++coord)
    ASSERT_DOUBLE_EQ (mapped_points[coord], expected[coord]);
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry_lagrange, LagrangeCmesh,
                          testing::Combine (testing::Range (T8_ECLASS_LINE, T8_ECLASS_TRIANGLE),
                                            testing::Range (1, MAX_POLYNOMIAL_DEGREE + 1)));

// TODO: use name_generator (http://google.github.io/googletest/reference/testing.html#INSTANTIATE_TEST_SUITE_P)
// to create test names like this: (EclassName_degree{degree}), e.g triangle_degree2, quad_degree3

TEST (foo, bar)
{
  t8_errorf ("\n-------------------\nFoo\n-------------------\n");
  std::vector<double> sp { 0.1, 1 };
  std::vector<double> nodes { 0, 0, 0, 0, 1, 0 };
  Boundary1D boundary (1);
  boundary.setNodes (nodes);
  auto pts = boundary.sample (sp);
  for (auto &i : pts)
    std::cout << i << ", ";
}
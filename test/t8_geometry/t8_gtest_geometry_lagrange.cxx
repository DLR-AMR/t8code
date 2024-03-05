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

#include <array>
#include <cstdlib>
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

#define MAX_POLYNOMIAL_DEGREE 2

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
allclose (const U &a, const V &b, double tol = 1e-8)
{
  size_t size_a = a.size ();
  size_t size_b = b.size ();
  if (size_a != size_b)
    return false;
  for (size_t i = 0; i < size_a; ++i)
    if (fabs (a[i] - b[i]) > tol)
      return false;
  return true;
}

/**
 * Flatten a vector of vector into a single vector.
 * 
 * \tparam T   Template parameter of a vector.
 * \param vec  Nested vector to be flattened.
 * \return     Flattened vector.
 */
template <typename T>
std::vector<T>
flatten (const std::vector<std::vector<T>> &vec)
{
  std::vector<T> flattened;
  for (auto const &v : vec) {
    flattened.insert (flattened.end (), v.begin (), v.end ());
  }
  return flattened;
}

t8_forest_t
create_uniform_forest (t8_cmesh_t cmesh, uint level)
{
  t8_forest_t forest;
  forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0, sc_MPI_COMM_WORLD);
  return forest;
};

/**
 * A single coarse mesh cell with Lagrange geometry.
 * 
 * This class is essentially a wrapper around a cmesh.
 * By having a single element instead of a mesh, understanding and
 * testing is made easier. Several topological utilities are provided,
 * some specific to the to the Lagrange geometry, some valid for all
 * the geometries in t8code.
 */
class LagrangeElement {
 public:
  /**
   * Construct a new LagrangeElement object.
   * 
   * \param eclass  Element class (line, quad, etc.)
   * \param degree  Polynomial degree (1, 2, ...)
   * \param nodes   x,y,z coordinates of the nodes, adhering to the numbering
   *                convention.
   */
  LagrangeElement (t8_eclass_t eclass, uint degree, std::vector<double> &nodes)
    : eclass (eclass), degree (degree), nodes (nodes)
  {
    // TODO: Check if the number of nodes corresponds to the element type and degree.
    // if (nodes.size () != parametric_nodes.size ())
    //   SC_ABORTF ("Provide the 3 coordinates of the nodes.\n");
    /* Create a cmesh with a single element */
    int dim = t8_eclass_to_dimension[eclass];
    t8_geometry_c *geometry = new t8_geometry_lagrange (dim);  // TODO: do it with a smart pointer or using RAII
    t8_cmesh_init (&cmesh);
    t8_cmesh_set_attribute (cmesh, 0, t8_get_package_id (), T8_CMESH_LAGRANGE_POLY_DEGREE, &degree, sizeof (int), 1);
    t8_cmesh_register_geometry (cmesh, geometry);
    t8_cmesh_set_tree_class (cmesh, 0, eclass);
    t8_cmesh_set_tree_vertices (cmesh, 0, nodes.data (), nodes.size ());
    t8_cmesh_commit (cmesh, sc_MPI_COMM_WORLD);
  }

  /**
   * Destroy the LagrangeElement object.
   * 
   * The cmesh wrapped by this class is also destroyed.
   * 
   */
  ~LagrangeElement ()
  {
    t8_cmesh_destroy (&cmesh);
  };

  /**
   * Get the type of the element.
   * 
   * \return  Element class of the element.
   */
  t8_eclass_t
  getType () const
  {
    return eclass;
  }

  /**
   * Element classes of the faces of this element.
   * 
   * Sandro: Not specific to the Lagrange geometry, may go into t8code base?
   * 
   * \return  Element classes of the faces, enumerated according to the face
   * ordering conventions of t8code.
   */
  std::vector<t8_eclass_t>
  faceClasses () const
  {
    std::vector<t8_eclass_t> face_classes;
    switch (eclass) {
    case T8_ECLASS_LINE:
      face_classes = std::vector<t8_eclass_t> (2, T8_ECLASS_VERTEX);
      break;
    case T8_ECLASS_QUAD:
      face_classes = std::vector<t8_eclass_t> (4, T8_ECLASS_LINE);
      break;
    case T8_ECLASS_TRIANGLE:
      face_classes = std::vector<t8_eclass_t> (3, T8_ECLASS_LINE);
      break;
    case T8_ECLASS_HEX:
      face_classes = std::vector<t8_eclass_t> (6, T8_ECLASS_QUAD);
      break;
    case T8_ECLASS_TET:
      face_classes = std::vector<t8_eclass_t> (4, T8_ECLASS_TRIANGLE);
      break;
    case T8_ECLASS_PRISM:
      face_classes = { T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE };
      break;
    case T8_ECLASS_PYRAMID:
      face_classes = { T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_TRIANGLE, T8_ECLASS_QUAD };
      break;
    default:
      SC_ABORTF ("Invalid element class %d.\n", eclass);
    }
    return face_classes;
  }

  /**
   * Coordinates of the specified node.
   * 
   * \param node  Node label. Node numbering starts at 0.
   * \return      x,y,z coordinates of the node.
   */
  std::vector<double>
  getNodeCoords (uint node) const
  {
    const double *v = t8_cmesh_get_tree_vertices (cmesh, 0);
    return std::vector<double> (v + 3 * node, v + 3 * node + 3);
  }

  /**
   * Coordinates of the specified nodes.
   * 
   * \param nodes  Node labels. Node numbering starts at 0.
   * \return       x,y,z coordinates of the nodes.
   */
  std::vector<std::vector<double>>
  getNodeCoords (std::vector<uint> &nodes) const
  {
    const double *v = t8_cmesh_get_tree_vertices (cmesh, 0);
    size_t n_node = nodes.size ();
    std::vector<std::vector<double>> node_coords (n_node);
    for (size_t i = 0; i < n_node; ++i) {
      uint i_node = nodes[i];
      node_coords[i] = std::vector<double> (v + 3 * i_node, v + 3 * i_node + 3);
    }
    return node_coords;
  }

  /**
   * Node labels on the faces of the element.
   * 
   * Sandro: What about making this function part of t8_geometry_lagrange?
   * 
   * \return  Node labels on each face of the element.
   */
  std::vector<std::vector<uint>>
  getFaceNodes () const
  {
    std::ostringstream invalid_degree;
    invalid_degree << "Invalid degree " << degree << ".\n";
    std::vector<std::vector<uint>> face_nodes;
    switch (eclass) {
    case T8_ECLASS_LINE:
      if (degree == 1)
        face_nodes = { { 0 }, { 1 } };
      else if (degree == 2)
        face_nodes = { { 0 }, { 1 }, { 2 } };
      else
        SC_ABORTF (invalid_degree.str ().c_str ());
      break;
    case T8_ECLASS_TRIANGLE:
      if (degree == 1)
        face_nodes = { { 1, 2 }, { 2, 0 }, { 0, 1 } };
      else if (degree == 2)
        face_nodes = { { 1, 2, 3 }, { 2, 0, 4 }, { 0, 1, 5 } };
      else
        SC_ABORTF (invalid_degree.str ().c_str ());
      break;
    case T8_ECLASS_QUAD:
      if (degree == 1) {
        face_nodes = { { 2, 0 }, { 1, 3 }, { 0, 1 }, { 3, 2 } };
      }
      else if (degree == 2)
        face_nodes = { { 2, 0, 4 }, { 1, 3, 5 }, { 0, 1, 6 }, { 3, 2, 7 } };
      else
        SC_ABORTF (invalid_degree.str ().c_str ());
      break;
    case T8_ECLASS_HEX:
      if (degree == 1) {
        face_nodes = { { 2, 0, 6, 4 }, { 1, 3, 5, 7 }, { 0, 1, 4, 5 }, { 3, 2, 7, 6 }, { 2, 3, 0, 1 }, { 4, 5, 6, 7 } };
      }
      else if (degree == 2) {
        face_nodes = { { 2, 0, 6, 4, 8, 9, 10, 11, 12 },   { 1, 3, 5, 7, 13, 14, 15, 16, 17 },
                       { 0, 1, 4, 5, 9, 13, 18, 19, 20 },  { 3, 2, 7, 6, 14, 8, 21, 22, 23 },
                       { 2, 3, 0, 1, 10, 15, 21, 18, 24 }, { 4, 5, 6, 7, 11, 16, 19, 22, 25 } };
      }
      else
        SC_ABORTF (invalid_degree.str ().c_str ());
      break;
    default:
      SC_ABORTF ("Invalid element class %d.\n", eclass);
    }
    return face_nodes;
  }

  /**
   * Decompose the element into its faces.
   * 
   * LagrangeElement()-s of codimension 1 are created. The original element
   * is not modified and the decomposition is not recursive.
   * The faces can further be decomposed by calling this method on them.
   * 
   * \return  Lagrange elements from the faces.
   * 
   * Sandro: If you find this decomposition useful, I could move it to
   * `t8_geometry_lagrange.hxx/cxx`.
   */
  std::vector<LagrangeElement>
  decompose () const
  {
    /* Get the node numbers of the faces */
    std::vector<t8_eclass_t> fc = faceClasses ();
    std::vector<std::vector<uint>> fn = getFaceNodes ();
    /* Create a new Lagrange element from each face */
    std::vector<LagrangeElement> faces;
    const uint n_face = t8_eclass_num_faces[eclass];
    faces.reserve (n_face);
    for (size_t i_face = 0; i_face < n_face; ++i_face) {
      auto nc = flatten<double> (getNodeCoords (fn[i_face]));
      faces.emplace_back (fc[i_face], degree, nc);
    }
    return faces;
  }

  /**
   * Physical coordinates of a point given in the reference domain.
   * 
   * \param point  Parametric coordinates of the point to be mapped.
   *               For 2D elements, only the first two coordinates are
   *               considered. For the 1D line element, only the first.
   * \return       Coordinates in the physical space.
   */
  std::array<double, 3>
  evaluate (const std::array<double, 3> &ref_point) const
  {
    std::array<double, 3> mapped;
    t8_geometry_evaluate (cmesh, 0, ref_point.data (), 1, mapped.data ());
    return mapped;
  }

  /**
   * Sample random points in the reference domain.
   * 
   * \param n_point  Number of points to generate.
   * \return         Coordinates of the points, given in x,y,z.
   */
  std::vector<std::array<double, 3>>
  sample (uint n_point) const
  {
    std::srand (time (NULL));
    std::vector<std::array<double, 3>> points (n_point);
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
   * Map this element on the face of a higher-dimensional element.
   * \verbatim
    Example to map a line element onto the left face of a quad element.

      Code: mapOnFace (T8_ECLASS_QUAD, 0, std::vector<double> { 0.6 })

      o : vertices of the elements
      + : point to map
                                y ^
                   x              |   face 3
    o------+--o  -->             2             3
          0.6                     o ----<---- o          faces (=edges in 2D)
                                  |           |          with CCW orientation
                                  |           |
                          face 0  v           ^  face 1
                                  + (0, 0.4)  |
                                  |           |   x
                                  o ---->---- o  -->
                                 0             1
                                      face 2
     \endverbatim
   * Sandro: If you find this mapping feauture useful, I could move it to
   * `t8_geometry_lagrange.hxx/cxx` and make it a public method (or a function).
   * 
   * \param eclass   Element class of the element onto which we map.
   *                 d-dimensional elements can only be mapped to the faces of
   *                 d+1-dimensional elements.
   * \param face_id  Face ID of the element onto which we map.
   *                 Faces in an element are numbered according to the t8code
   *                 conventions. The selected face must have the same type as
   *                 the element from which we map. For instance, an element of
   *                 type T8_ECLASS_QUAD can only be mapped onto face 4 of a
   *                 T8_ECLASS_PYRAMID, since the other faces of a pyramidal
   *                 element are triangles.
   * \param coord    x,y,z coordinates of the point in this element.
   * \return         x,y,z coordinates of the mapped point.
   */
  std::array<double, 3>
  mapOnFace (t8_eclass map_onto, const uint face_id, const std::array<double, 3> &coord) const
  {
    /* Error messages for input validation */
    std::ostringstream unsupported_element;
    unsupported_element << "Mapping from a " << t8_eclass_to_string[eclass] << " element is not supported yet.\n";
    std::ostringstream unsupported_target_element;
    unsupported_target_element << "Mapping onto a " << t8_eclass_to_string[map_onto]
                               << " element is not supported yet.\n";
    std::ostringstream too_many_faces;
    too_many_faces << "A " << t8_eclass_to_string[map_onto] << " element has " << t8_eclass_num_faces[map_onto]
                   << " faces only.\n";
    std::ostringstream non_matching_elem_classes;
    non_matching_elem_classes << "Face " << face_id << " of a " << t8_eclass_to_string[map_onto]
                              << " element is not of type " << t8_eclass_to_string[eclass] << ".\n";
    if (face_id > t8_eclass_num_faces[map_onto] - 1)
      SC_ABORT (too_many_faces.str ().c_str ());

    /* Actual mapping, case by case */
    std::array<double, 3> mapped_coord;
    double xi = coord[0];
    double eta = coord[1];
    double zeta = coord[2];
    switch (eclass) {
    case T8_ECLASS_LINE:
      if (map_onto == T8_ECLASS_TRIANGLE) {
        if (face_id == 0)
          mapped_coord = { 1, xi, 0 };
        else if (face_id == 1)
          mapped_coord = { 1 - xi, 1 - xi, 0 };
        else if (face_id == 2)
          mapped_coord = { xi, 0, 0 };
      }
      else if (map_onto == T8_ECLASS_QUAD) {
        if (face_id == 0)
          mapped_coord = { 0, 1 - xi, 0 };
        else if (face_id == 1)
          mapped_coord = { 1, xi, 0 };
        else if (face_id == 2)
          mapped_coord = { xi, 0, 0 };
        else if (face_id == 3)
          mapped_coord = { 1 - xi, 1, 0 };
      }
      else
        SC_ABORT (unsupported_target_element.str ().c_str ());
      break;
    case T8_ECLASS_QUAD:
      if (map_onto == T8_ECLASS_HEX) {
        if (face_id == 0)
          mapped_coord = { 0, 1 - xi, eta };
        else if (face_id == 1)
          mapped_coord = { 1, xi, eta };
        else if (face_id == 2)
          mapped_coord = { xi, 0, eta };
        else if (face_id == 3)
          mapped_coord = { 1 - xi, 1, eta };
        else if (face_id == 4)
          mapped_coord = { xi, 1 - eta, 0 };
        else if (face_id == 5)
          mapped_coord = { xi, eta, 1 };
      }
      else
        SC_ABORT (unsupported_target_element.str ().c_str ());
      break;
    default:
      SC_ABORT (unsupported_element.str ().c_str ());
    }
    return mapped_coord;
  }

  /**
   * Save the geometry into a VTK file.
   * 
   * \remark Note that the geometry is exported as a linear cell for now,
   * so the exported result is accurate for \a degree 1 only.
   * 
   */
  void
  write () const
  {
    /* A cmesh cannot be exported, only a forest.
       So we create one tree element per coarse mesh element */
    t8_forest_t forest = create_uniform_forest (cmesh, 0);
    std::ostringstream filename;
    filename << "Lagrange" << t8_eclass_to_string[eclass] << "Degree" << degree;
    t8_forest_write_vtk (forest, filename.str ().c_str ());
    /* Clean up the dummy forest, making sure that the cmesh is not destroyed */
    t8_cmesh_ref (cmesh);
    t8_forest_unref (&forest);
  }

  /**
   * Create a sample LagrangeElement.
   * 
   * The goal of this factory method is to quickly instantiate a
   * LagrangeElement for the purpose of testing.
   * 
   * \param eclass  Element class of the element.
   * \param degree  Polynomial degree.
   * \return        LagrangeElement.
   */
  static LagrangeElement
  create_sample_element (t8_eclass_t eclass, int degree)
  {
    std::ostringstream invalid_degree;
    invalid_degree << "Degree " << degree << " is not yet supported for " << t8_eclass_to_string[eclass]
                   << " elements.\n";
    std::vector<double> vertices;
    double perturb_amplitude = 1 / 5.0 * (1 / (double) degree);
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
        vertices = perturb (vertices, perturb_amplitude);
        break;
      case 2:
        vertices = { 0,   0, 0,   1,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   1, 1,   0,   1,   0,   1, 1,
                     1,   1, 1,   0,   1,   0.5, 0,   0,   0.5, 0,   0.5, 0,   0,   0.5, 1, 0,   0.5, 0.5, 1,   0, 0.5,
                     1,   1, 0.5, 1,   0.5, 0,   1,   0.5, 1,   1,   0.5, 0.5, 0.5, 0,   0, 0.5, 0,   1,   0.5, 0, 0.5,
                     0.5, 1, 0,   0.5, 1,   1,   0.5, 1,   0.5, 0.5, 0.5, 0,   0.5, 0.5, 1, 0.5, 0.5, 0.5 };
        vertices = perturb (vertices, perturb_amplitude);
        break;
      default:
        SC_ABORT (invalid_degree.str ().c_str ());
      }
      break;
    default:
      SC_ABORTF ("Not implemented for %s elements.\n", t8_eclass_to_string[eclass]);
    }
    return LagrangeElement (eclass, degree, vertices);
  }

 private:
  /** Lagrange elements have the same element class as the linear ones. */
  t8_eclass_t eclass;
  /** Polynomial degree of the geometrical mapping. */
  const uint degree;
  /** Points in the physical space, which span the geometry of the element. */
  const std::vector<double> nodes;
  /** Coarse mesh, wrapped by this class. */
  t8_cmesh_t cmesh;
  /** Number of nodes in the Lagrange element of a given class and degree */
  static constexpr uint lagrange_nodes[T8_ECLASS_COUNT][2] = { { 1, 1 }, { 2, 3 }, { 4, 9 }, { 3, 6 }, { 8, 27 } };
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
  LagrangeElement lag = LagrangeElement::create_sample_element (eclass, degree);
  lag.write ();
  t8_errorf ("\n-------------------\nEclass: %s, degree: %d\n-------------------\n", t8_eclass_to_string[eclass],
             degree);
  /* Compare the mappings on each boundary face of the Lagrange element */
  std::vector<LagrangeElement> faces = lag.decompose ();
  uint i_face = 0;
  for (const auto &face : faces) {
    auto points_on_face = face.sample (5);
    for (const auto &point : points_on_face) {
      auto mapped1 = face.evaluate (point);
      auto mapped2 = lag.evaluate (face.mapOnFace (lag.getType (), i_face, point));
      ASSERT_TRUE (allclose (mapped1, mapped2));
    }
    ++i_face;
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_geometry_lagrange, LagrangeCmesh,
                          testing::Combine (testing::Range (T8_ECLASS_QUAD, T8_ECLASS_TET),
                                            testing::Range (1, MAX_POLYNOMIAL_DEGREE + 1)),
                          [] (const testing::TestParamInfo<LagrangeCmesh::ParamType> &info) {
                            std::ostringstream test_name;
                            test_name << t8_eclass_to_string[std::get<0> (info.param)] << "_degree"
                                      << std::get<1> (info.param);
                            return test_name.str ();
                          });

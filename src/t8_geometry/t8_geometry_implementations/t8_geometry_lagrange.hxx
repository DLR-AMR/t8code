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

/** \file t8_geometry_lagrange.hxx
 * This geometry implements curved geometries obtained by Lagrange interpolation.
 * The interpolation is carried out with the classical finite element basis functions.
 * Therefore, it is the generalization of the linear geometry (t8_geometry_linear).
 */

#ifndef T8_GEOMETRY_LAGRANGE_HXX
#define T8_GEOMETRY_LAGRANGE_HXX

#include <array>
#include <string>
#include <vector>

#include <t8.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>

/** The maximum polynomial degree currently implemented. */
#define T8_GEOMETRY_MAX_POLYNOMIAL_DEGREE 2

/**
 * Mapping with Lagrange basis functions
 *
 * The enumeration of the nodal basis functions depends on the node numbering
 * scheme. While this is a convention, we came up with the following rules to
 * guide us when constructing new mappings:
 * 1. Be compatible with lower degree elements. It means that when we read the
 * first \a n nodes, it gives a valid lower-degree element. This rule is made
 * attainable by the fact that the nodes spanning the Lagrange basis are
 * equidistant. As an example, consider a degree four segment element whose
 * five nodes in increasing coordinate direction are numbered as 0-3-2-4-1.
 * Reading the first two nodes gives us a valid linear segment, while reading
 * two first three nodes provides a valid quadratic segment. In general, we
 * need to find for an element of degree \a d which nodes coincide with the
 * nodes of elements with degree lower than \a d. Those nodes must be numbered
 * first.
 * 2. Faces are oriented outwards the volume (3D) or the plane (2D). This is
 * just a convention to ensure that the node numbering be consistent.
 * 3. The mapping by the Lagrange geometry falls back to the linear geometry
 * for degree one. It means that the starting point for the node number
 * assignment is the numbering defined in the \a t8_element.cxx file.
 * 4. The node numbering is performed in increasing spatial dimension, which
 * results in a hierarchical construction of elements. The node numbers are
 * assigned based on increasing face IDs, then increasing edge IDs, then
 * increasing vertex IDs. For volume elements, the volume nodes are assigned
 * lastly. This hierarchical construction satisfies rule 1 at the mesh level.
 *
 * You can verify these rules by checking the documentation of the
 * \a t8_geom_xxx_basis methods of this class
 * (e.g. t8_geometry_lagrange::t8_geom_t6_basis).
 *
 */
struct t8_geometry_lagrange: public t8_geometry_with_vertices
{
 public:
  /**
   * Constructor of the Lagrange geometry with a given dimension. The geometry
   * is compatible with all tree types and uses as many vertices as the number of Lagrange
   * basis functions used for the mapping.
   * The vertices are saved via the \ref t8_cmesh_set_tree_vertices function.
   */
  t8_geometry_lagrange ();

  virtual ~t8_geometry_lagrange ();

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_LAGRANGE;
  };

  /**
   * Maps points from the reference space to the physical space \f$ \mathbb{R}^3 \f$.
   *
   * For linear elements, it gives the same result as \ref t8_geom_compute_linear_geometry.
   *
   * The mapping is performed via Lagrange interpolation according to
   *
   * \f$ \mathbf{x}(\vec{\xi}) = \sum\limits_{i=1}^{N_{\mathrm{vertex}}} \psi_i(\vec{\xi}) \mathbf{x}_i \f$
   *
   * where \f$ \vec{\xi} \f$ is the point in the reference space to be mapped, \f$ \mathbf{x} \f$ is the mapped point we search,
   * \f$ \psi_i(\vec{\xi}) \f$ are the basis functions associated with the vertices, and \f$ \mathbf{x}_i \f$ are the
   * vertices of the current tree in the physical space.
   * The basis functions are specific to the tree type, see e.g. t8_geom_t6_basis .
   * The vertices of the current tree were set with \ref t8_cmesh_set_tree_vertices.
   *
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_points entries, specifying points in the reference space.
   * \param [in]  num_points  Number of points to map. Currently, only one point is supported.
   * \param [out] out_coords  Coordinates of the mapped points in physical space of \a ref_coords. The length is \a num_points * 3.
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_points,
                    double *out_coords) const;

  /**
   * Compute the Jacobian of the \a t8_geom_evaluate map at a point in the reference space.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_points entries, specifying points in the reference space.
   * \param [in]  num_points  Number of points to map.
   * \param [out] jacobian    The Jacobian at \a ref_coords. Array of size \a num_points x dimension x 3. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the Jacobian  (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   */
  virtual void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_points,
                             double *jacobian) const;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the polynomial degree of the
   * Lagrange basis into an internal buffer (as is done in the Lagrange geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * This geometry supports lines, triangles, quadrilaterals and hexahedra up to degree 2.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const;

 private:
  /**
   * Evaluates the basis functions of the current tree type at a point.
   * \param [in]  ref_point  Array of tree dimension entries, specifying the point in the reference space.
   */
  inline std::vector<double>
  t8_geom_compute_basis (const double *ref_point) const;

  /**
   * Basis functions of a 2-node segment.
   * \verbatim
      x --------- x
     0             1
     \endverbatim
   * \param [in] ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_s2_basis (const double *ref_point) const;

  /**
   * Basis functions of a 3-node segment.
   * \verbatim
      x ----x---- x
     0      2      1
     \endverbatim
   * \param [in] ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_s3_basis (const double *ref_point) const;

  /**
   * Basis functions of a 3-node triangle element.
   * \verbatim
                   2
                  x
                / |
              /   |
            /     |
          /       |
        /         |
      x --------- x
     0             1
     \endverbatim
   * \param [in] ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_t3_basis (const double *ref_point) const;

  /**
   * Basis functions of a 6-node triangle element.
   * \verbatim
                   2
                  x
                / |
              /   |
          4 x     x 3
          /       |
        /         |
      x --- x --- x
     0      5      1
     \endverbatim
   * \param [in] ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_t6_basis (const double *ref_point) const;

  /**
   * Basis functions of a 4-node quadrilateral element.
   * \verbatim
     2             3
      x --------- x
      |           |
      |           |
      |           |
      |           |
      |           |
      x --------- x
     0             1
     \endverbatim
   * \param [in] ref_point Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_q4_basis (const double *ref_point) const;

  /**
   * Basis functions of a 9-node quadrilateral element.
   * \verbatim
     2      7      3
      x ----x---- x
      |           |
      |     8     |
     4x     x     x5
      |           |
      |           |
      x ----x---- x
     0      6      1
     \endverbatim
   * \param [in] ref_point Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_q9_basis (const double *ref_point) const;

  /**
   * Basis functions of an 8-node hexahedron element.
   * \verbatim
     z
      ^   6           7
      |  x ----------x
        /|          /|
     4 / |       5 / |
      x --------- x  |
      |  |        |  |
      |  |2       |  |3
      |  x--------|--x
      | /         | /
      |/          |/     x
      x --------- x    -->
     0             1
     \endverbatim
   * \param [in] ref_point Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_h8_basis (const double *ref_point) const;

  /**
   * Basis functions of a 27-node hexahedron element.
   * \verbatim
          z = 0                     z = 0.5                   z = 1
     y                         y                         y
      ^                         ^                         ^
      |                         |                         |
     2     21     3            8     23      14          6     22      7
      x ----x---- x             x ----x---- x             x ----x---- x
      |           |             |           |             |           |
      |    24     |             |    26     |             |    25     |
    10x     x     x15         12x     x     x17         11x     x     x16
      |           |             |           |             |           |
      |           |    x        |           |    x        |           |    x
      x ----x---- x  -->        x ----x---- x  -->        x ----x---- x  -->
     0     18     1            9     20      13          4     19      5
     \endverbatim
   * \param ref_point [in] Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  inline std::vector<double>
  t8_geom_h27_basis (const double *ref_point) const;

  /** Polynomial degree of the interpolation. */
  const int *degree;
};

/**
 * Flatten a vector of vector into a single vector.
 *
 * \tparam T   Template parameter of a vector.
 * \param [in] vec Nested vector to be flattened.
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

/**
 * A single coarse mesh cell with Lagrange geometry.
 *
 * This class is essentially a wrapper around a cmesh.
 * By having a single element instead of a mesh, understanding and
 * testing is made easier. Several topological utilities are provided,
 * some specific to the Lagrange geometry, some valid for all the
 * geometries in t8code.
 */
struct t8_lagrange_element
{
 public:
  /**
   * Construct a new t8_lagrange_element object.
   *
   * \param [in] eclass  Element class (line, quad, etc.)
   * \param [in] degree  Polynomial degree (1, 2, ...)
   * \param [in] nodes   x,y,z coordinates of the nodes, adhering to the numbering
   *                     convention.
   */
  t8_lagrange_element (t8_eclass_t eclass, uint32_t degree, std::vector<double> &nodes);

  /**
   * Destroy the t8_lagrange_element object.
   *
   * The cmesh wrapped by this class is also destroyed.
   *
   */
  ~t8_lagrange_element ()
  {
    t8_cmesh_destroy (&cmesh);
  };

  /**
   * Get the type of the element.
   *
   * \return  Element class of the element.
   */
  t8_eclass_t
  get_type () const;

  /**
   * Element classes of the faces of this element.
   *
   * \return  Element classes of the faces, enumerated according to the face
   * ordering conventions of t8code.
   */
  std::vector<t8_eclass_t>
  face_classes () const;

  /**
   * Coordinates of the specified node.
   *
   * \param [in] node  Node label. Node numbering starts at 0.
   * \return           x,y,z coordinates of the node.
   */
  std::vector<double>
  get_node_coords (uint32_t node) const;

  /**
   * Coordinates of the specified nodes.
   *
   * \param [in] nodes  Node labels. Node numbering starts at 0.
   * \return            x,y,z coordinates of the nodes.
   */
  std::vector<std::vector<double>>
  get_node_coords (std::vector<uint32_t> &nodes) const;

  /**
   * Node labels on the faces of the element.
   *
   * \return  Node labels on each face of the element.
   */
  std::vector<std::vector<uint32_t>>
  get_face_nodes () const;

  /**
   * Decompose the element into its faces.
   *
   * t8_lagrange_element()-s of codimension 1 are created. The original element
   * is not modified and the decomposition is not recursive.
   * The faces can further be decomposed by calling this method on them.
   *
   * \return  Lagrange elements from the faces.
   */
  std::vector<t8_lagrange_element>
  decompose () const;

  /**
   * Physical coordinates of a point given in the reference domain.
   *
   * \param [in] ref_point  Parametric coordinates of the point to be mapped.
   *                        For 2D elements, only the first two coordinates are
   *                        considered. For the 1D line element, only the first.
   * \return                Coordinates in the physical space.
   */
  std::array<double, T8_ECLASS_MAX_DIM>
  evaluate (const std::array<double, T8_ECLASS_MAX_DIM> &ref_point) const;

  /**
   * Sample random points in the reference domain.
   *
   * \param [in] n_point  Number of points to generate.
   * \return              Coordinates of the points, given in x,y,z.
   */
  std::vector<std::array<double, T8_ECLASS_MAX_DIM>>
  sample (uint32_t n_point) const;

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
   *
   * \param [in] eclass   Element class of the element onto which we map.
   *                      d-dimensional elements can only be mapped to the faces of
   *                      d+1-dimensional elements.
   * \param [in] face_id  Face ID of the element onto which we map.
   *                      Faces in an element are numbered according to the t8code
   *                      conventions. The selected face must have the same type as
   *                      the element from which we map. For instance, an element of
   *                      type \ref T8_ECLASS_QUAD can only be mapped onto face 4 of a
   *                      \ref T8_ECLASS_PYRAMID, since the other faces of a pyramidal
   *                      element are triangles.
   * \param [in] coord    x,y,z coordinates of the point in this element.
   * \return              x,y,z coordinates of the mapped point.
   */
  std::array<double, T8_ECLASS_MAX_DIM>
  map_on_face (t8_eclass eclass, const int face_id, const std::array<double, T8_ECLASS_MAX_DIM> &coord) const;

  /**
   * Save the geometry into a VTK file.
   *
   * \remark Note that the geometry is exported as a linear cell for now,
   * so the exported result is accurate for \a degree 1 only.
   *
   */
  void
  write () const;

 private:
  /** Lagrange elements have the same element class as the linear ones. */
  t8_eclass_t eclass;
  /** Polynomial degree of the geometrical mapping. */
  const uint32_t degree;
  /** Points in the physical space, which span the geometry of the element. */
  const std::vector<double> nodes;
  /** Coarse mesh, wrapped by this class. */
  t8_cmesh_t cmesh;
  /** Number of nodes in the Lagrange element of a given class and degree */
  static constexpr uint32_t lagrange_nodes[T8_ECLASS_COUNT][2] = { { 1, 1 }, { 2, 3 }, { 4, 9 }, { 3, 6 }, { 8, 27 } };

  /** Inbuilt function to create a uniform forest.
   * \param [in] cmesh The cmesh to use.
   * \param [in] level The level of the uniform forest.
   * \return           The forest.
   */
  t8_forest_t
  create_uniform_forest (t8_cmesh_t cmesh, uint32_t level) const;
};

#endif /* !T8_GEOMETRY_LAGRANGE_HXX */

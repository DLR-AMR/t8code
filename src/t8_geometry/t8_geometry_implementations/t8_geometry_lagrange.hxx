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

/** \file t8_geometry_lagrange.hxx
 * This geometry implements curved geometries obtained by Lagrange interpolation.
 * The interpolation is carried out with the classical finite element basis functions.
 * Therefore, it is the generalization of the linear geometry (t8_geometry_linear).
 */

#ifndef T8_GEOMETRY_LAGRANGE_HXX
#define T8_GEOMETRY_LAGRANGE_HXX

#include <string>
#include <vector>

#include <t8.h>
#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>

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
 * assignment is the numbering defined in the \a t8_element.c file.
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
   * \param [in] dim  0 <= \a dim <= 3. Element dimension in the parametric space.
   * E.g. \a dim = 2 for a \ref T8_ECLASS_QUAD element.
   */
  t8_geometry_lagrange (int dim);

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets dimension and name to invalid values. */
  t8_geometry_lagrange (): t8_geometry_with_vertices ()
  {
  }

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
   * For linear elements, it gives the same result as \ref t8_geom_compute_linear_geometry.
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_points entries, specifying points in the reference space.
   * \param [in]  num_points  Number of points to map.
   * \param [out] out_coords  Coordinates of the mapped points in physical space of \a ref_coords. The length is \a num_points * 3.
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_points,
                    double *out_coords) const;

  /**
   * Compute the Jacobian of the \a t8_geom_evaluate map at a point in the reference space.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_points entries, specifying points in the reference space.
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

 private:
  /**
   * Evaluates the basis functions of the current tree type at a point.
   * \param [in]  ref_point  Array of \a dimension entries, specifying the point in the reference space.
   */
  const std::vector<double>
  t8_geom_compute_basis (const double *ref_point) const;

  /**
   * Map a point from the reference space to the physical space.
   * The mapping is performed via Lagrange interpolation according to
   * 
   * \f$ \mathbf{x}(\vec{\xi}) = \sum\limits_{i=1}^{N_{\mathrm{vertex}}} \psi_i(\vec{\xi}) \mathbf{x}_i \f$
   * 
   * where \f$ \vec{\xi} \f$ is the point in the reference space to be mapped, \f$ \mathbf{x} \f$ is the mapped point we search,
   * \f$ \psi_i(\vec{\xi}) \f$ are the basis functions associated with the vertices, and \f$ \mathbf{x}_i \f$ are the
   * vertices of the current tree in the physical space.
   * The basis functions are specific to the tree type, see e.g. \ref t8_geom_t6_basis.
   * The vertices of the current tree were set with \ref t8_cmesh_set_tree_vertices.
   * \param [in]  ref_point     Array of \a dimension entries, specifying the coordinates of \f$ \vec{\xi} \f$.
   * \param [out] mapped_point  Array of 3 entries, specifying the coordinates of \f$ \mathbf{x} \f$.
   */
  void
  t8_geom_map (const double *ref_point, double *mapped_point) const;

  /**
   * Basis functions of a 2-node segment.
   * \verbatim
      x --------- x
     0             1
     \endverbatim
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
  t8_geom_s2_basis (const double *ref_point) const;

  /**
   * Basis functions of a 3-node segment.
   * \verbatim
      x ----x---- x
     0      2      1
     \endverbatim
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
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
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
  t8_geom_t3_basis (const double *ref_point) const;

  /**
   * Basis functions of a 6-node triangle element.
   * \verbatim
                   2
                  x
                / |
              /   |
          5 x     x 4
          /       |
        /         |
      x --- x --- x
     0      3      1
     \endverbatim
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
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
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
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
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
  t8_geom_q9_basis (const double *ref_point) const;

  /**
   * Basis functions of an 8-node hexahedron element.
   * \verbatim
     y
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
   * \param ref_point  Point in the reference space.
   * \return  Basis functions evaluated at the reference point.
   */
  const std::vector<double>
  t8_geom_h8_basis (const double *ref_point) const;

  /** Polynomial degree of the interpolation. */
  const int *degree;
};

#endif /* !T8_GEOMETRY_LAGRANGE_HXX! */

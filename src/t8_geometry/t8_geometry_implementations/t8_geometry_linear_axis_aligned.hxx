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

/** \file t8_geometry_linear_axis_aligned.hxx
 * Definition of an axis-aligned geometry. It maps from
 * \f$ [0,1]^\mathrm{dim} \f$ to \f$ \mathbb{R}^3 \f$, using two vertices. 
 */

#ifndef T8_GEOMETRY_LINEAR_AXIS_ALIGNED_HXX
#define T8_GEOMETRY_LINEAR_AXIS_ALIGNED_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>

struct t8_geometry_linear_axis_aligned: public t8_geometry_with_vertices
{
 public:
  /**
   * Constructor of the linear, axis-aligned geometry with a given dimension.
   * The geometry is only viable for the line/quad/hex tree types and uses two
   * vertices (min and max coords) per tree. The vertices are saved via
   * the \ref t8_cmesh_set_tree_vertices function. Sets the dimension and the
   * name to "t8_geom_linear_axis_aligned_{dim}"
   * \param [in] dim  0 <= \a dimension <= 3. The dimension.
   */
  t8_geometry_linear_axis_aligned (int dim);

  /* Base constructor with no arguments. We need this since it
   * is called from derived class constructors.
   * Sets dimension and name to invalid values. */

  t8_geometry_linear_axis_aligned (): t8_geometry_with_vertices ()
  {
  }

  /** The destructor. 
   */
  virtual ~t8_geometry_linear_axis_aligned ();

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const
  {
    return T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED;
  };

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   */
  virtual void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \a num_coords x dimension x 3. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian  (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   */
  virtual void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *jacobian) const;

  /**
   * \param[in] forest            The forest of the element.
   * \param[in] ltreeid           The local tree id of the element's tree
   * \param[in] element           The element
   * \param[in] points            points to check
   * \param[in] num_points        Number of points to check
   * \param[in, out] is_inside    Array to fill with flags whether the point is inside or not
   * \param[in] tolerance         Tolerance of the inside-check
   */
  virtual void
  t8_geom_point_batch_inside_element (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                      const double *points, const int num_points, int *is_inside,
                                      const double tolerance) const;

  /**
   * Check if the currently active tree has a negative volume.
   * \return                True if the currently loaded tree has a negative volume.
   */
  virtual bool
  t8_geom_tree_negative_volume () const;

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * Only line, quad and hex elements are supported by this geometry.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const override
  {
    if (active_tree_class == T8_ECLASS_LINE || active_tree_class == T8_ECLASS_QUAD
        || active_tree_class == T8_ECLASS_HEX) {
      t8_productionf ("Axis-aligned geometry is not compatible with tree type %s\n It is only compatible with line, "
                      "quad and hex elements.\n",
                      t8_eclass_to_string[active_tree_class]);
      return true;
    }
    return false;
  }
};

#endif /* !T8_GEOMETRY_LINEAR_AXIS_ALIGNED_HXX */

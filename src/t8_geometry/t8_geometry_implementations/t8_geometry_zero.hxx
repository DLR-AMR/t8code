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

/** \file t8_geometry_zero.hxx
 * This geometry implements the t8_geometry interface and maps
 * all points to 0.
 * This geometry serves for testing and debugging purposes.
 */

#ifndef T8_GEOMETRY_ZERO_HXX
#define T8_GEOMETRY_ZERO_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>
#include <t8_geometry/t8_geometry_helpers.h>

struct t8_geometry_zero: public t8_geometry
{
 public:
  /**
   * Constructor of the zero geometry. The geometry
   * is viable with all tree types. This geometry maps all points to zero and
   * is meant for debugging purposes.
   * Sets the dimension and the name to "t8_geom_zero"
   */
  t8_geometry_zero ();

  /**
   * Check if  the currently active tree has a negative volume
   * \return                True (non-zero) if the currently loaded tree has a negative volume. 0 otherwise.  
   */
  bool
  t8_geom_tree_negative_volume () const override
  {
    return 0;
  };

  /** The destructor. 
   * Clears the allocated memory.
   */
  virtual ~t8_geometry_zero ();

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const override
  {
    return T8_GEOMETRY_TYPE_ZERO;
  };

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
   * \note All entries in out_coords will be set to 0.
   */
  void
  t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                    double *out_coords) const override;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \a num_coords x dimension x 3. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian  (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   * \note All entries in \a jacobian will be set to zero.
   */
  void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *jacobian) const override;

  /**
   * \param[in] forest            The forest of the element.
   * \param[in] ltreeid           The local tree id of the element's tree
   * \param[in] element           The element
   * \param[in] points            points to check
   * \param[in] num_points        Number of points to check
   * \param[in, out] is_inside    Array to fill with flags whether the point is inside or not
   * \param[in] tolerance         Tolerance of the inside-check
   */
  void
  t8_geom_point_batch_inside_element (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                      const double *points, const int num_points, int *is_inside,
                                      const double tolerance) const override
  {
    const double zeros[3] = { 0 };
    for (int i_point = 0; i_point < num_points; ++i_point) {
      const int offset = i_point * T8_ECLASS_MAX_DIM;
      is_inside[i_point] = t8_vertex_point_inside (zeros, points + offset, tolerance);
    }
  }

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated.
   * In this implementation, we do nothing, since we do not need any tree data.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  inline void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid) override;

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * This geometry supports all element types, hence it returns true.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const override
  {
    return true;
  }
};

#endif /* T8_GEOMETRY_ZERO_HXX */

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

/** \file t8_geometry_analytic.hxx
 * This geometry implements analytic geometries. It provides an interface to
 * define custom functions for evaluation, Jacobians and the loading of the
 * tree data.
 */

#ifndef T8_GEOMETRY_ANALYTIC_HXX
#define T8_GEOMETRY_ANALYTIC_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_analytic.h>

struct t8_geometry_analytic: public t8_geometry
{
 public:
  /**
   * Constructor of the analytic geometry. The geometry
   * is viable with all tree types and uses a user-provided analytic and
   * jacobian function. The actual mappings are done by these functions.
   * \param [in] name       The name to give this geometry.
   * \param [in] analytical The analytical function to use for this geometry.
   * \param [in] jacobian   The jacobian of \a analytical.
   * \param [in] load_tree_data The function that is used to load a tree's data.
   * \param [in] tree_negative_volume_in The function that is used to compute if a trees volume is negative.
   * \param [in] tree_compatible_in The function that is used to check if a tree is compatible with the geometry.
   */
  t8_geometry_analytic (std::string name, t8_geom_analytic_fn analytical, t8_geom_analytic_jacobian_fn jacobian,
                        t8_geom_load_tree_data_fn load_tree_data,
                        t8_geom_tree_negative_volume_fn tree_negative_volume_in,
                        t8_geom_tree_compatible_fn tree_compatible_in, const void *user_data);

  /**
   * Constructor of the analytic geometry for testing purposes.
   * \param [in] name       The name to give this geometry.
   */
  t8_geometry_analytic (std::string name);

  /** The destructor. 
   */
  virtual ~t8_geometry_analytic ()
  {
    /* Nothing to do */
  }

  /**
   * Get the type of this geometry.
   * \return The type.
   */
  inline t8_geometry_type_t
  t8_geom_get_type () const override
  {
    return T8_GEOMETRY_TYPE_ANALYTIC;
  };

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of tree dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of \f$ \mathrm{dim} \f$ to map.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
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
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \f$ \mathrm{dim} \cdot 3 \f$ x \a num_coords. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   * \note The jacobian will be
   *            (1)              (1 0)             (1 0 0)
   * dim 1: J = (0)   dim 2: J = (0 1)  dim 3: J = (0 1 0)
   *            (0)              (0 0)             (0 0 1)
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
  t8_geom_point_batch_inside_element ([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] t8_locidx_t ltreeid, [[maybe_unused]] const t8_element_t *element,
    [[maybe_unused]] const double *points, [[maybe_unused]] const int num_points, [[maybe_unused]] int *is_inside,
    [[maybe_unused]] const double tolerance) const override
  {
    SC_ABORTF ("Function not yet implemented");
  }

  /**
   * Check if the currently active tree has a negative volume
   * \return                True (non-zero) if the currently loaded tree has a negative volume. 0 otherwise.  
   */
  bool
  t8_geom_tree_negative_volume () const override;

  /**
   * Check for compatibility of the currently loaded tree with the geometry.
   * If the geometry has limitations these can be checked here.
   * This includes for example if only specific tree types or dimensions are supported.
   * If all trees are supported, this function should return true.
   * \return                True if the geometry is compatible with the tree.
   */
  bool
  t8_geom_check_tree_compatibility () const override;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the vertex coordinates of the 
   * tree into an internal buffer (as is done in the linear geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid) override;

  inline const void *
  t8_geom_analytic_get_user_data ()
  {
    return user_data;
  }

 private:
  t8_geom_analytic_fn analytical_function; /**< The given analytical function. */

  t8_geom_analytic_jacobian_fn jacobian; /**< Its jacobian. */

  t8_geom_load_tree_data_fn load_tree_data; /**< The function to load the tree data. */

  t8_geom_tree_negative_volume_fn tree_negative_volume; /**< The function to check for negative volumes. */

  t8_geom_tree_compatible_fn tree_compatible; /**< The function to check if a tree is compatible. */

  const void *tree_data; /** Tree data pointer that can be set in \a load_tree_data and 
                             is passed onto \a analytical_function and \a jacobian. */

  const void *user_data; /** Additional user data pointer that can be set in constructor
                             and modified via \ref t8_geom_analytic_get_user_data. */
};

#endif /* !T8_GEOMETRY_ANALYTICAL_HXX */

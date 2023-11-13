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
 * TODO: Add description
 */

#ifndef T8_GEOMETRY_ANALYTIC_HXX
#define T8_GEOMETRY_ANALYTIC_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_with_vertices.hxx>
#include <t8_geometry/t8_geometry_with_vertices.h>

/**
 * Definition of an analytic geometry function.
 * This function maps reference coordinates to physical
 * coordinates.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of dimension x \a num_coords many entries, specifying a point in \f$ [0,1]^\mathrm{dim} \f$.
 * \param [in]  num_coords  
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords. The length is \a num_coords * 3.
 * \param [in]  tree_data   The data of the current tree as loaded by a \ref t8_geom_load_tree_data_fn.
 * \param [in]  user_data   The user data pointer stored in the geometry.
 */
typedef void (*t8_geom_analytic_fn) (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                     const size_t num_coords, double *out_coords, const void *tree_data,
                                     const void *user_data);

/**
 * Definition for the jacobian of an analytic geometry function.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
 * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
 * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \f$ \mathrm{dim} \cdot 3 \f$ x \a num_coords. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
 *                          correspond to the \f$ i \f$-th column of the jacobian (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
 * \param [in]  tree_data   The data of the current tree as loaded by a \ref t8_geom_load_tree_data_fn.
 * \param [in]  user_data   The user data pointer stored in the geometry.
 */
typedef void (*t8_geom_analytic_jacobian_fn) (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                              const size_t num_coords, double *jacobian, const void *tree_data,
                                              const void *user_data);

/* TODO: Document. */
typedef void (*t8_geom_load_tree_data_fn) (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const void **tree_data);

struct t8_geometry_analytic: public t8_geometry
{
 public:
  /**
   * Constructor of the analytic geometry with a given dimension. The geometry
   * is viable with all tree types and uses a user-provided analytic and
   * jacobian function. The actual mappings are done by these functions.
   * \param [in] dim        The dimension of this geometry.
   * \param [in] name       The name to give this geometry.
   * \param [in] analytical The analytical function to use for this geometry.
   * \param [in] jacobian   The jacobian of \a analytical.
   * \param [in] load_tree_data The function that is used to load a tree's data.
   */
  t8_geometry_analytic (int dim, const char *name, t8_geom_analytic_fn analytical,
                        t8_geom_analytic_jacobian_fn jacobian, t8_geom_load_tree_data_fn load_tree_data,
                        const void *user_data);

  /** The destructor. 
   * Clears the allocated memory.
   */
  virtual ~t8_geometry_analytic ()
  {
    /* Nothing to do */
  }

  /**
   * Maps points in the reference space \f$ [0,1]^\mathrm{dim} \to \mathbb{R}^3 \f$.
   * \param [in]  cmesh       The cmesh in which the point lies.
   * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension x \a num_coords many entries, specifying points in \f$ [0,1]^\mathrm{dim} \f$.
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
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
   * \param [in]  num_coords  Amount of points of /f$ \mathrm{dim} /f$ to map.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size \f$ \mathrm{dim} \cdot 3 \f$ x \a num_coords. Indices \f$ 3 \cdot i\f$ , \f$ 3 \cdot i+1 \f$ , \f$ 3 \cdot i+2 \f$
   *                          correspond to the \f$ i \f$-th column of the jacobian (Entry \f$ 3 \cdot i + j \f$ is \f$ \frac{\partial f_j}{\partial x_i} \f$).
   * \note The jacobian will be
   *            (1)              (1 0)             (1 0 0)
   * dim 1: J = (0)   dim 2: J = (0 1)  dim 3: J = (0 1 0)
   *            (0)              (0 0)             (0 0 1)
   */
  virtual void
  t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                             double *jacobian) const;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated. You can use it for example to load the vertex coordinates of the 
   * tree into an internal buffer (as is done in the linear geometry).
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual void
  t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid);

  inline const void *
  t8_geom_analytic_get_user_data ()
  {
    return user_data;
  }

 private:
  t8_geom_analytic_fn analytical_function; /**< The given analytical function. */

  t8_geom_analytic_jacobian_fn jacobian; /**< Its jacobian. */

  t8_geom_load_tree_data_fn load_tree_data; /**< The function to load the tree data. */

  const void *tree_data; /** Tree data pointer that can be set in \a load_tree_data and 
                                           is passed onto \a analytical_function and \a jacobian. */

  const void *user_data; /** Additional user data pointer that can be set in constructor
                                         * and modified via \ref t8_geom_analytic_get_user_data. */
};

/* TODO: Document */
void
t8_geom_load_tree_data_vertices (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const void **vertices_out);

#endif /* !T8_GEOMETRY_ANALYTICAL_HXX! */

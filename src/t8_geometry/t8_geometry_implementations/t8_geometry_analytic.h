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

/** \file t8_geometry_analytic.h
 * This header provides the C interface to create an analytical geometry.
 */

#ifndef T8_GEOMETRY_ANALYTIC_H
#define T8_GEOMETRY_ANALYTIC_H

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

/**
 * Definition for the load tree data function.
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree (of the cmesh) in which the reference point is.
 * \param [in]  tree_data   The data of the trees.
 */
typedef void (*t8_geom_load_tree_data_fn) (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const void **tree_data);

/**
 * Definition for the negative volume function.
 */
typedef bool (*t8_geom_tree_negative_volume_fn) ();

T8_EXTERN_C_BEGIN ();

/** Destroy a geometry analytic object.
 * \param [in,out] geom   A pointer to a geometry object. Set to NULL on output.
 */
void
t8_geometry_analytic_destroy (t8_geometry_c **geom);

/** Create a new analytical geometry.
 * \return          A pointer to an allocated geometry struct.
 */
t8_geometry_c *
t8_geometry_analytic_new (int dim, const char *name, t8_geom_analytic_fn analytical,
                          t8_geom_analytic_jacobian_fn jacobian, t8_geom_load_tree_data_fn load_tree_data,
                          t8_geom_tree_negative_volume_fn tree_negative_volume, const void *user_data);

/**
 * Load vertex data from given tree. 
 * \param [in]  cmesh       The cmesh.
 * \param [in]  gtreeid     The global tree id (in the cmesh).
 * \param [out] vertex_out  The load tree vertices.
 */
void
t8_geom_load_tree_data_vertices (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const void **user_data);

T8_EXTERN_C_END ();

#endif /* T8_GEOMETRY_ANALYTIC_H */

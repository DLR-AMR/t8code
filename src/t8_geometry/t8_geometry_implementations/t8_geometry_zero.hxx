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
 * This geometry servers for testing and debugging purposes.
 */

#ifndef T8_GEOMETRY_ZERO_HXX
#define T8_GEOMETRY_ZERO_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>

struct t8_geometry_zero:public t8_geometry
{
public:

  /* Basic constructor that sets the dimension and the name
   * to "t8_geom_zero_{dimension}" */
  t8_geometry_zero (int dimension);

  /** The destructor. 
   * Clears the allocated memory.
   */
                      virtual ~ t8_geometry_zero ();

  /**
   * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   * \note All entries in out_coords will be set to 0.
   */
  virtual void        t8_geom_evaluate (t8_cmesh_t cmesh,
                                        t8_gloidx_t gtreeid,
                                        const double *ref_coords,
                                        double out_coords[3]) const;

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  cmesh      The cmesh in which the point lies.
   * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   * \note All entries in \a jacobian will be set to zero.
   */
  virtual void        t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                                t8_gloidx_t gtreeid,
                                                const double *ref_coords,
                                                double *jacobian) const;

  /** Update a possible internal data buffer for per tree data.
   * This function is called before the first coordinates in a new tree are
   * evaluated.
   * In this implementation, we do nothing, since we do not need any tree data.
   * \param [in]  cmesh      The cmesh.
   * \param [in]  gtreeid    The global tree.
   */
  virtual inline void t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid);
};

#endif /* !T8_GEOMETRY_ZERO_HXX! */

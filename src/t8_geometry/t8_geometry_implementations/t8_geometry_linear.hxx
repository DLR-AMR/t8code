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

/** \file t8_geometry_linear.hxx
 * TODO: Add description
 */

#ifndef T8_GEOMETRY_LINEAR_HXX
#define T8_GEOMETRY_LINEAR_HXX

#include <t8.h>
#include <t8_geometry/t8_geometry_base.hxx>

struct t8_geometry_linear:public t8_geometry
{
public:

  /* Basic constructor that sets the dimension and the name
   * to "t8_geom_linear_{dimension}" */
  t8_geometry_linear (int dimension);

  /** The destructor. 
   * Clears the allocated memory.
   */
                      virtual ~ t8_geometry_linear ();

  /**
   * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
   * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
   */
  virtual void        t8_geom_evaluate (t8_gloidx_t ltree_id,
                                        const double *ref_coords,
                                        double out_coords[3]);

  /**
   * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
   * \param [in]  ltree_id    The local tree (of the cmesh) in which the reference point is.
   * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
   * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
   *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
   */
  virtual void        t8_geom_evalute_jacobian (t8_gloidx_t ltree_id,
                                                const double *ref_coords,
                                                double *jacobian);

  virtual inline void t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid);

protected:
                      t8_gloidx_t active_tree;
  t8_eclass_t         active_tree_class;
  const double       *active_tree_vertices;
};

#endif /* !T8_GEOMETRY_LINEAR_HXX! */

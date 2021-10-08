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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>

t8_geometry_zero::t8_geometry_zero (int dim)
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t              num_chars = 100;
  char               *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_zero_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_zero::~t8_geometry_zero ()
{
  T8_FREE ((char *) name);
}

/**
 * Map a point in the reference space $$[0,1]^dimension$$ to $$\mathbb R^3$$
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] out_coords  The mapped coordinates in physical space of \a ref_coords.
 */
void
t8_geometry_zero::t8_geom_evaluate (t8_cmesh_t cmesh,
                                    t8_gloidx_t gtreeid,
                                    const double *ref_coords,
                                    double out_coords[3]) const
{
  /* Set the out_coords to 0 */
  memset (out_coords, 0, sizeof (double) * 3);
}

/**
 * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
 *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
 */
void
t8_geometry_zero::t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                            t8_gloidx_t gtreeid,
                                            const double
                                            *ref_coords,
                                            double *jacobian) const
{
  /* Set the jacobian to 0 */
  memset (jacobian, 0, sizeof (double) * 3 * dimension);
}

inline void
t8_geometry_zero::t8_geom_load_tree_data (t8_cmesh_t cmesh,
                                          t8_gloidx_t gtreeid)
{
  /* Do nothing. */
}

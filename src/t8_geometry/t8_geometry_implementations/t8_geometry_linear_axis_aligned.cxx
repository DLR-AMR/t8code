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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear_axis_aligned.h>
#include <t8_geometry/t8_geometry_helpers.h>

t8_geometry_linear_axis_aligned::t8_geometry_linear_axis_aligned (int dim):
t8_geometry_w_vertices (dim, "")
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t              num_chars = 100;
  char               *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_linear_axis_aligned_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_linear_axis_aligned::~t8_geometry_linear_axis_aligned ()
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
/* *INDENT-OFF* */
/* indent adds second const */
void
t8_geometry_linear_axis_aligned::t8_geom_evaluate (t8_cmesh_t cmesh,
                                      t8_gloidx_t gtreeid,
                                      const double *ref_coords,
                                      double out_coords[3]) const
/* *INDENT-ON* */
{
  t8_geom_compute_linear_axis_aligned_geometry (active_tree_class,
                                                active_tree_vertices,
                                                ref_coords, out_coords);
}

/**
 * Compute the jacobian of the \a t8_geom_evaluate map at a point in the reference space $$[0,1]^dimension$$.
 * \param [in]  cmesh      The cmesh in which the point lies.
 * \param [in]  gtreeid    The global tree (of the cmesh) in which the reference point is.
 * \param [in]  ref_coords  Array of \a dimension many entries, specifying a point in [0,1]^dimension.
 * \param [out] jacobian    The jacobian at \a ref_coords. Array of size dimension x 3. Indices 3*i, 3*i+1, 3*i+2
 *                          correspond to the i-th column of the jacobian (Entry 3*i + j is del f_j/del x_i).
 */
/* *INDENT-OFF* */
/* indent adds second const */
void
t8_geometry_linear_axis_aligned::t8_geom_evalute_jacobian (t8_cmesh_t cmesh,
                                              t8_gloidx_t gtreeid,
                                              const double
                                              *ref_coords,
                                              double *jacobian) const
/* *INDENT-ON* */
{
  SC_ABORT ("Not implemented.");
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear_axis_aligned.h.
 * Create a new geometry with given dimension. */
t8_geometry_c      *
t8_geometry_linear_axis_aligned_new (int dimension)
{
  t8_geometry_linear_axis_aligned *geom =
    new t8_geometry_linear_axis_aligned (dimension);
  return (t8_geometry_c *) geom;
}

void
t8_geometry_linear_axis_aligned_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () ==
             T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED);

  delete             *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();

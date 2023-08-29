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

t8_geometry_linear_axis_aligned::t8_geometry_linear_axis_aligned (int dim): t8_geometry_with_vertices (dim, "")
{
  T8_ASSERT (0 <= dim && dim <= 3);
  size_t num_chars = 100;
  char *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_linear_axis_aligned_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_linear_axis_aligned::~t8_geometry_linear_axis_aligned ()
{
  T8_FREE ((char *) name);
}

void
t8_geometry_linear_axis_aligned::t8_geom_evaluate (t8_cmesh_t cmesh,
                                                   t8_gloidx_t gtreeid,
                                                   const double *ref_coords,
                                                   const size_t num_coords,
                                                   double out_coords[3]) const
{
  t8_geom_compute_linear_axis_aligned_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords,
                                                out_coords);
}

void
t8_geometry_linear_axis_aligned::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh,
                                                            t8_gloidx_t
                                                            gtreeid,
                                                            const double
                                                            *ref_coords,
                                                            const size_t
                                                            num_coords,
                                                            double *jacobian)
  const
{
  SC_ABORT ("Not implemented.");
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear_axis_aligned.h.
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_linear_axis_aligned_new (int dimension)
{
  t8_geometry_linear_axis_aligned *geom = new t8_geometry_linear_axis_aligned (dimension);
  return (t8_geometry_c *) geom;
}

void
t8_geometry_linear_axis_aligned_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT (t8_geom_is_linear_axis_aligned (*geom));

  delete *geom;
  *geom = NULL;
}

#if T8_ENABLE_DEBUG
int
t8_geom_is_linear_axis_aligned (const t8_geometry_c *geometry)
{
  /* Try to dynamic cast the geometry into linear, axis-aligned geometry. 
   * This is only successful if geometry pointed to a 
   * t8_geometry_linear_axis_aligned.
   * If successful, then is_linear_geom will be true.
   */
  const int is_linear_axis_aligned_geom = (dynamic_cast<const t8_geometry_linear_axis_aligned *> (geometry) != NULL);

  return is_linear_axis_aligned_geom;
}
#endif

T8_EXTERN_C_END ();

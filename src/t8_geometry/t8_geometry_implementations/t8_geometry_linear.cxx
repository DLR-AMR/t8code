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

#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#include <t8_geometry/t8_geometry_helpers.h>

t8_geometry_linear::t8_geometry_linear (int dim)
  : t8_geometry_with_vertices (dim, "t8_geom_linear_" + std::to_string (dim))
{
}

t8_geometry_linear::~t8_geometry_linear ()
{
}

void
t8_geometry_linear::t8_geom_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                      const size_t num_coords, double *out_coords) const
{
  t8_geom_compute_linear_geometry (active_tree_class, active_tree_vertices, ref_coords, num_coords, out_coords);
}

void
t8_geometry_linear::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords,
                                               const size_t num_coords, double *jacobian) const
{
  SC_ABORT ("Not implemented.");
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_linear.h.
 * Create a new geometry with given dimension. */
t8_geometry_c *
t8_geometry_linear_new (int dimension)
{
  t8_geometry_linear *geom = new t8_geometry_linear (dimension);
  return (t8_geometry_c *) geom;
}

void
t8_geometry_linear_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_LINEAR);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();

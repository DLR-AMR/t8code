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

/** \file t8_geometry_zero.cxx
 * Implements functions declared in \ref t8_geometry_zero.hxx.
 */

#include <t8_geometry/t8_geometry_implementations/t8_geometry_zero.hxx>

t8_geometry_zero::t8_geometry_zero (): t8_geometry ("t8_geom_zero")
{
}

t8_geometry_zero::~t8_geometry_zero ()
{
}

void
t8_geometry_zero::t8_geom_evaluate ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                                    [[maybe_unused]] const double *ref_coords, const size_t num_coords,
                                    double *out_coords) const
{
  /* Set the out_coords to 0 */
  for (size_t coord = 0; coord < num_coords; coord++) {
    out_coords[0 + num_coords * T8_ECLASS_MAX_DIM] = 0;
    out_coords[1 + num_coords * T8_ECLASS_MAX_DIM] = 0;
    out_coords[2 + num_coords * T8_ECLASS_MAX_DIM] = 0;
  }
}

void
t8_geometry_zero::t8_geom_evaluate_jacobian ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid,
                                             [[maybe_unused]] const double *ref_coords, const size_t num_coords,
                                             double *jacobian) const
{
  /* Set the jacobian to 0 */
  const int tree_dim = t8_eclass_to_dimension[active_tree_class];
  memset (jacobian, 0, sizeof (double) * 3 * tree_dim * num_coords);
}

inline void
t8_geometry_zero::t8_geom_load_tree_data ([[maybe_unused]] t8_cmesh_t cmesh, [[maybe_unused]] t8_gloidx_t gtreeid)
{
  /* Do nothing. */
}

T8_EXTERN_C_BEGIN ();

/* Satisfy the C interface from t8_geometry_zero.h.
 * Create a new geometry. */
t8_geometry_c *
t8_geometry_zero_new ()
{
  t8_geometry_zero *geom = new t8_geometry_zero ();
  return (t8_geometry_c *) geom;
}

void
t8_geometry_zero_destroy (t8_geometry_c **geom)
{
  T8_ASSERT (geom != NULL);
  T8_ASSERT ((*geom)->t8_geom_get_type () == T8_GEOMETRY_TYPE_ZERO);

  delete *geom;
  *geom = NULL;
}

T8_EXTERN_C_END ();

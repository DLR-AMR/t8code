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
  size_t num_chars = 100;
  char *name_tmp = T8_ALLOC (char, num_chars);

  snprintf (name_tmp, num_chars, "t8_geom_zero_%i", dim);
  name = name_tmp;
  dimension = dim;
}

t8_geometry_zero::~t8_geometry_zero ()
{
  T8_FREE ((char *) name);
}

void
t8_geometry_zero::t8_geom_evaluate (t8_cmesh_t cmesh,
                                    t8_gloidx_t gtreeid,
                                    const double *ref_coords,
                                    const size_t num_coords,
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
t8_geometry_zero::t8_geom_evaluate_jacobian (t8_cmesh_t cmesh,
                                             t8_gloidx_t gtreeid,
                                             const double
                                             *ref_coords,
                                             const size_t num_coords,
                                             double *jacobian) const
{
  /* Set the jacobian to 0 */
  memset (jacobian, 0, sizeof (double) * 3 * dimension * num_coords);
}

inline void
t8_geometry_zero::t8_geom_load_tree_data (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  /* Do nothing. */
}

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

#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry.h>
#include <t8_geometry/t8_geometry_handler.hxx>

void
t8_geometry_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *out_coords)
{
  double start_wtime = 0; /* Used for profiling. */
  /* The geometries do not expect the in- and output vector to be the same */
  T8_ASSERT (ref_coords != out_coords);

  if (cmesh->profile != NULL) {
    /* Measure the runtime of geometry evaluation.
     * We accumulate the runtime over all calls. */
    start_wtime = sc_MPI_Wtime ();
  }

  if (cmesh->geometry_handler == NULL) {
    SC_ABORT ("Error: Trying to evaluate non-existing geometry.\n");
  }

  /* Evaluate the geometry. */
  cmesh->geometry_handler->evaluate_tree_geometry (cmesh, gtreeid, ref_coords, num_coords, out_coords);

  if (cmesh->profile != NULL) {
    /* If profiling is enabled, add the runtime to the profiling
     * variable. */
    cmesh->profile->geometry_evaluate_runtime += sc_MPI_Wtime () - start_wtime;
    cmesh->profile->geometry_evaluate_num_calls++;
  }
}

void
t8_geometry_jacobian (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *jacobian)
{
  /* Evaluate the jacobian. */
  cmesh->geometry_handler->evaluate_tree_geometry_jacobian (cmesh, gtreeid, ref_coords, num_coords, jacobian);
}

t8_geometry_type_t
t8_geometry_get_type (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  if (cmesh->geometry_handler == nullptr) {
    return T8_GEOMETRY_TYPE_INVALID;
  }
  /* Return the type. */
  return cmesh->geometry_handler->get_tree_geometry_type (cmesh, gtreeid);
}

int
t8_geometry_tree_negative_volume (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid)
{
  return cmesh->geometry_handler->tree_negative_volume (cmesh, gtreeid);
}

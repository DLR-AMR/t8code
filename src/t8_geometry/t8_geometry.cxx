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

#include <t8_cmesh/t8_cmesh_types.h>
#include <t8_geometry/t8_geometry.h>

void
t8_geometry_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *out_coords)
{
  double start_wtime = 0; /* Used for profiling. */
  /* The cmesh must be committed */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* The geometries do not expect the in- and output vector to be the same */
  T8_ASSERT (ref_coords != out_coords);
  /* Get the geometry handler of the cmesh. */
  t8_geometry_handler &geom_handler = cmesh->geometry_handler;

  if (cmesh->profile != NULL) {
    /* Measure the runtime of geometry evaluation.
     * We accumulate the runtime over all calls. */
    start_wtime = sc_MPI_Wtime ();
  }
  /* Detect whether we call this function for the first time in a row for 
   * this tree and if so update the active tree and geometry. */
  geom_handler.update_tree (cmesh, gtreeid);

  /* Evaluate the geometry. */
  geom_handler.evaluate_active_geometry (cmesh, ref_coords, num_coords, out_coords);

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
  /* The cmesh must be committed */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* Get the geometry handler of the cmesh of the forest. */
  t8_geometry_handler &geom_handler = cmesh->geometry_handler;

  /* Detect whether we call this function for the first time in a row for 
   * this tree and if so update the active tree and geometry. */
  geom_handler.update_tree (cmesh, gtreeid);

  /* Evaluate the jacobian. */
  geom_handler.evaluate_active_geometry_jacobian (cmesh, ref_coords, num_coords, jacobian);
}

t8_geometry_type_t
t8_geometry_get_type (t8_cmesh_t cmesh, t8_gloidx_t gtreeid)
{
  /* The cmesh must be committed */
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  /* Get the geometry handler of the cmesh of the forest. */
  t8_geometry_handler &geom_handler = cmesh->geometry_handler;

  /* Detect whether we call this function for the first time in a row for 
   * this tree and if so update the active tree and geometry. */
  geom_handler.update_tree (cmesh, gtreeid);

  /* Return the type. */
  return geom_handler.get_active_geometry_type ();
}

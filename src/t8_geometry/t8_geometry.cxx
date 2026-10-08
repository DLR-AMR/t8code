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
#include <t8_geometry/t8_geometry.hxx>
#include <t8_geometry/t8_geometry_c_interface.h>
#include <t8_geometry/t8_geometry_handler.hxx>
#include <algorithm>
#include <vector>

void
t8_geometry_evaluate (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, std::span<const t8_3D_vec> ref_coords,
                      std::span<t8_3D_vec> out_coords)
{
  T8_ASSERT (ref_coords.size () == out_coords.size ());
  const size_t num_coords = ref_coords.size ();
  if (num_coords == 0) {
    return;
  }
  if (num_coords == 1) {
    out_coords[0] = t8_geometry_evaluate (cmesh, gtreeid, ref_coords[0]);
    return;
  }
  /** Number of entries per point that the geometries expect in their reference coordinates input, max (d, 1). */
  const size_t ref_dim = static_cast<size_t> (std::max (t8_cmesh_get_dimension (cmesh), 1));
  /** The reference coordinates as consecutive points with \a ref_dim entries each, as expected by the geometries. */
  std::vector<double> packed_ref_coords (num_coords * ref_dim);
  for (size_t icoord = 0; icoord < num_coords; ++icoord) {
    std::copy_n (ref_coords[icoord].begin (), ref_dim, packed_ref_coords.begin () + icoord * ref_dim);
  }
  /** The domain coordinates as consecutive points with 3 entries each, as returned by the geometries. */
  std::vector<double> packed_out_coords (3 * num_coords);
  t8_geometry_evaluate (cmesh, gtreeid, packed_ref_coords.data (), num_coords, packed_out_coords.data ());
  for (size_t icoord = 0; icoord < num_coords; ++icoord) {
    std::copy_n (packed_out_coords.begin () + 3 * icoord, 3, out_coords[icoord].begin ());
  }
}

t8_3D_vec
t8_geometry_evaluate (const t8_cmesh_t cmesh, const t8_gloidx_t gtreeid, const t8_3D_vec &ref_coords)
{
  /* A single point does not need to be repacked, since the geometry only reads its first d components. */
  t8_3D_vec out_coords;
  t8_geometry_evaluate (cmesh, gtreeid, ref_coords.data (), 1, out_coords.data ());
  return out_coords;
}

T8_EXTERN_C_BEGIN ();

void
t8_geometry_evaluate (t8_cmesh_t cmesh, t8_gloidx_t gtreeid, const double *ref_coords, const size_t num_coords,
                      double *out_coords)
{
  double start_wtime = 0; /* Used for profiling. */
  /* The geometries do not expect the in- and output vector to be the same */
  T8_ASSERT (ref_coords != out_coords);

  if (cmesh->profile != nullptr) {
    /* Measure the runtime of geometry evaluation.
     * We accumulate the runtime over all calls. */
    start_wtime = sc_MPI_Wtime ();
  }

  if (cmesh->geometry_handler == nullptr) {
    SC_ABORT ("Error: Trying to evaluate non-existing geometry.\n");
  }

  /* Evaluate the geometry. */
  cmesh->geometry_handler->evaluate_tree_geometry (cmesh, gtreeid, ref_coords, num_coords, out_coords);

  if (cmesh->profile != nullptr) {
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

T8_EXTERN_C_END ();

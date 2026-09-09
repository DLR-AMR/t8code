/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_vertex_conn_helpers.cxx
 * Routines that derive vertex connectivity of an uncommitted cmesh from its
 * trees' corner coordinates.
 */

#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_helpers.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_helpers.hxx>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_connectivity.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

#include <array>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cstring>

std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>>
t8_cmesh_vertex_conn_helpers_global_ids_by_coordinates (const t8_stash_t stash,
                                                        const std::vector<t8_eclass_t> &eclasses)
{
  const t8_gloidx_t ntrees = stash->classes.elem_count;
  const std::vector<std::vector<t8_3D_vec>> vertices = t8_stash_extract_vertices (stash);
  SC_CHECK_ABORT (vertices.empty (), "t8_cmesh_vertex_conn_helpers_global_ids_by_coordinates: cmesh has no tree "
                                     "vertex coordinates set; cannot derive global vertex ids from coordinates.");

  const t8_cmesh_coord_match_context ctx = t8_cmesh_compute_coord_match_context (ntrees, eclasses.data (), vertices);

  std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS> unset;
  unset.fill (-1);
  std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>> global_ids (ntrees, unset);
  std::unordered_multimap<unsigned long, std::pair<int, int>> vertex_bins;
  t8_gloidx_t next_id = 0;

  for (int itree = 0; itree < ntrees; itree++) {
    const int nverts = t8_eclass_num_vertices[eclasses[itree]];

    for (int ivert = 0; ivert < nverts; ivert++) {
      const t8_3D_vec &vertex = vertices[itree][ivert];
      const unsigned long hash = t8_cmesh_hash_coord (ctx, vertex);

      t8_gloidx_t matched_id = -1;
      auto range = vertex_bins.equal_range (hash);
      for (auto it = range.first; it != range.second; ++it) {
        const int cand_itree = it->second.first;
        const int cand_ivert = it->second.second;
        const t8_3D_vec &cand_vertex = vertices[cand_itree][cand_ivert];

        if (t8_cmesh_coords_match (ctx, vertex, cand_vertex)) {
          matched_id = global_ids[cand_itree][cand_ivert];
          break;
        }
      }

      global_ids[itree][ivert] = (matched_id > -1) ? matched_id : next_id++;
      vertex_bins.insert (std::make_pair (hash, std::make_pair (itree, ivert)));
    }
  }

  return global_ids;
}

void
t8_cmesh_vertex_conn_set_vertices_by_coordinates (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!t8_cmesh_is_committed (cmesh, 0));

  const t8_stash_t stash = cmesh->stash;
  const t8_gloidx_t ntrees = stash->classes.elem_count;

  if (ntrees == 0) {
    /* Nothing to derive ids for. */
    return;
  }

  const std::vector<t8_eclass_t> eclasses = t8_stash_extract_eclasses (stash);

  const std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>> global_ids
    = t8_cmesh_vertex_conn_helpers_global_ids_by_coordinates (stash, eclasses);

  t8_cmesh_enable_vertex_conn (cmesh);
  for (int itree = 0; itree < ntrees; itree++) {
    const int nverts = t8_eclass_num_vertices[eclasses[itree]];
    t8_cmesh_set_global_vertices_of_tree (cmesh, itree, global_ids[itree].data (), nverts);
  }
}

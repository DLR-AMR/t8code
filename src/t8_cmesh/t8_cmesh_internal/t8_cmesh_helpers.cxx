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

/** \file t8_cmesh_helpers.cxx
 *
 * Collection of cmesh helper routines.
 */

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_helpers.hxx>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>
#include <t8_cmesh/t8_cmesh_vertex_connectivity/t8_cmesh_vertex_conn_helpers.hxx>

#include <array>
#include <cmath>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <utility>

t8_cmesh_coord_match_context
t8_cmesh_compute_coord_match_context (const t8_gloidx_t ntrees, const t8_eclass_t *eclasses,
                                      const std::vector<std::vector<t8_3D_vec>> &vertices)
{
  /* Compute minimum and maximum of the cmesh domain. */
  double min_coord = vertices[0][0][0];
  double max_coord = vertices[0][0][0];

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];
    const int nverts = t8_eclass_num_vertices[eclass];
    const int edim = t8_eclass_to_dimension[eclass];

    for (int ivert = 0; ivert < nverts; ivert++) {
      const t8_3D_vec &vertex = vertices[itree][ivert];

      for (int icoord = 0; icoord < edim; icoord++) {
        const double coord = vertex[icoord];

        if (coord < min_coord) {
          min_coord = coord;
        }

        if (coord > max_coord) {
          max_coord = coord;
        }
      }
    }
  }

  /* Scaled tolerance to decide if two doubles are equal. */
  const double tolerance = 10.0 * T8_PRECISION_EPS * std::abs (max_coord - min_coord);

  /* `num_bins` should be more than enough for (almost) all cases.
   * I.e., 2^P4EST_QMAXLEVEL =~ 1.073e9.
   */
  const double num_bins = 1e9; /* Number of bins. */
  const double inverse_bin_size = num_bins / (max_coord - min_coord);

  return { min_coord, max_coord, tolerance, inverse_bin_size };
}

unsigned long
t8_cmesh_hash_coord (const t8_cmesh_coord_match_context &ctx, const t8_3D_vec &vertex)
{
  /* Compute the hash key. The idea is to convert the rescaled vertex coordinates to
   * long integers and add them up. We apply a bit of seeding by also adding `icoord`. */
  unsigned long hash = 0;

  for (int icoord = 0; icoord < T8_ECLASS_MAX_DIM; icoord++) {
    const double rescaled = (vertex[icoord] - ctx.min_coord) * ctx.inverse_bin_size;

    /* Simple hash function. */
    hash = hash + icoord + static_cast<unsigned long> (rescaled + 0.5);
  }

  return hash;
}

bool
t8_cmesh_coords_match (const t8_cmesh_coord_match_context &ctx, const t8_3D_vec &vertex_a, const t8_3D_vec &vertex_b)
{
  /* All T8_ECLASS_MAX_DIM coordinate components must agree within `ctx.tolerance`. */
  return t8_eq (vertex_a, vertex_b, ctx.tolerance);
}

/** Fill \a global_ids (resized to one fixed-size corner array per tree) with the global vertex
 * ids of every tree currently in \a stash that has any set.
 * \return True if at least one tree's global vertex ids were found in \a stash. */
static bool
t8_cmesh_helpers_extract_global_ids (const t8_stash_t stash, const t8_gloidx_t ntrees,
                                     std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>> &global_ids)
{
  std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS> unset;
  unset.fill (-1);
  global_ids.assign (ntrees, unset);
  bool found_global_ids = false;

  for (size_t iattribute = 0; iattribute < stash->attributes.elem_count; iattribute++) {
    const t8_stash_attribute_struct_t *entry
      = (const t8_stash_attribute_struct_t *) t8_sc_array_index_locidx (&stash->attributes, iattribute);
    if (entry->key == T8_CMESH_GLOBAL_VERTICES_ATTRIBUTE_KEY && entry->package_id == t8_get_package_id ()) {
      memcpy (global_ids[entry->id].data (), entry->attr_data, entry->attr_size);
      found_global_ids = true;
    }
  }

  return found_global_ids;
}

/** Find faces whose vertices all correspond (have the same global vertex id in \a global_ids) to
 * another face's vertices, and join them. Shared core of \ref t8_cmesh_join_by_vertices.
 */
static void
t8_cmesh_join_by_vertices_impl (const t8_cmesh_t cmesh, const t8_gloidx_t ntrees,
                                const std::vector<t8_eclass_t> &eclasses,
                                const std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>> &global_ids)
{
  /* Setup hash table `faces` mapping a hash key to a pair containing `(itree, iface)`. */
  std::unordered_multimap<unsigned long, std::pair<int, int>> faces;
  const std::vector<std::array<std::optional<std::tuple<t8_gloidx_t, int, int>>, T8_ECLASS_MAX_FACES>>
    face_already_joined = t8_stash_extract_joined_faces (cmesh->stash);

  for (int itree = 0; itree < ntrees; itree++) {
    const t8_eclass_t eclass = eclasses[itree];
    const int nfaces = t8_eclass_num_faces[eclass];

    for (int iface = 0; iface < nfaces; iface++) {
      if (face_already_joined[itree][iface].has_value ()) {
        /* This face already has an explicit join (e.g. a periodic boundary); leave it alone. */
        continue;
      }

      const int nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[eclass][iface]];

      /* Compute the hash key as the sum of the face's vertices' global ids. */
      unsigned long hash = 0;
      for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
        const int ivert = t8_face_vertex_to_tree_vertex[eclass][iface][iface_vert];
        hash += static_cast<unsigned long> (global_ids[itree][ivert]);
      }

      /* Loop over all pre-registered faces with the same hash. */
      auto range = faces.equal_range (hash);
      for (auto it = range.first; it != range.second; ++it) {
        const int neigh_itree = it->second.first;
        const int neigh_iface = it->second.second;

        const t8_eclass_t neigh_eclass = eclasses[neigh_itree];
        const int neigh_nface_verts = t8_eclass_num_vertices[t8_eclass_face_types[neigh_eclass][neigh_iface]];

        if (nface_verts != neigh_nface_verts) {
          continue;
        }

        /* The order of the encountered face vertices is needed for computing
         * the orientation later on. Prepare the array for that here. */
        int face_vert_order[T8_ECLASS_MAX_EDGES_2D];
        for (int i = 0; i < T8_ECLASS_MAX_EDGES_2D; i++) {
          face_vert_order[i] = -1;
        }

        int match_count = 0; /* This tracks the number of matching vertices. */
        for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
          const int ivert = t8_face_vertex_to_tree_vertex[eclass][iface][iface_vert];

          for (int neigh_iface_vert = 0; neigh_iface_vert < neigh_nface_verts; neigh_iface_vert++) {
            const int neigh_ivert = t8_face_vertex_to_tree_vertex[neigh_eclass][neigh_iface][neigh_iface_vert];

            if (global_ids[itree][ivert] == global_ids[neigh_itree][neigh_ivert]) {
              match_count++;
              /* Store the encountered face vertex order for later use. */
              face_vert_order[iface_vert] = neigh_iface_vert;
              continue;
            }
          }
        }

        /* If the number of matching face vertices is equal to the actual number of the face's
         * vertices we interpret this as a face-to-face connection between two elements. */
        if (match_count == nface_verts) {
          /* Compute the orientation of the face-to-face connection.
           * Face corner 0 of the face with the lower face direction connects
           * to a corner of the other face. The number of this corner is the
           * orientation code. */
          int orientation = -1;
          int smaller_bigger_face_condition = -1;

          const int compare = t8_eclass_compare (eclass, neigh_eclass);
          if (compare < 0) {
            /* This tree class is smaller than neigh. tree class. */
            smaller_bigger_face_condition = 1;
          }
          else if (compare > 0) {
            /* This tree class is bigger than neigh. tree class. */
            smaller_bigger_face_condition = 0;
          }
          else {
            /* This tree class is the same as the neigh. tree class.
               Then the face with the smaller face id is the smaller one. */
            smaller_bigger_face_condition = iface < neigh_iface;
          }

          if (smaller_bigger_face_condition) {
            orientation = face_vert_order[0];
          }
          else {
            for (int iface_vert = 0; iface_vert < nface_verts; iface_vert++) {
              if (0 == face_vert_order[iface_vert]) {
                orientation = iface_vert;
                break;
              }
            }
          }

          t8_cmesh_set_join (cmesh, itree, neigh_itree, iface, neigh_iface, orientation);

          break;
        }
      } /* Loop over faces with identical hash. */

      /* Register the current pair of `itree` and `iface` with given `hash` in the hash table. */
      faces.insert (std::make_pair (hash, std::make_pair (itree, iface)));
    } /* Loop over faces. */
  }   /* Loop over trees. */
}

void
t8_cmesh_join_by_vertices (const t8_cmesh_t cmesh)
{
  T8_ASSERT (t8_cmesh_is_initialized (cmesh));
  T8_ASSERT (!t8_cmesh_is_committed (cmesh, 0));

  const t8_stash_t stash = cmesh->stash;
  const t8_gloidx_t ntrees = stash->classes.elem_count;

  if (ntrees == 0) {
    /* Nothing to join. */
    return;
  }

  const std::vector<t8_eclass_t> eclasses = t8_stash_extract_eclasses (stash);

  /* Fast path: if global vertex ids are already present in the stash (either set manually by the
   * user or derived by t8_cmesh_vertex_conn_set_vertices_by_coordinates), reuse them to match
   * faces by exact id equality. Otherwise derive ids from coordinates on the fly, purely as an
   * internal speedup for this matching step: they are never persisted, so a caller that only
   * asked for automatic face joining (and not for vertex connectivity) sees no side effects. */
  std::vector<std::array<t8_gloidx_t, T8_ECLASS_MAX_CORNERS>> global_ids;
  if (!t8_cmesh_helpers_extract_global_ids (stash, ntrees, global_ids)) {
    global_ids = t8_cmesh_vertex_conn_helpers_global_ids_by_coordinates (stash, eclasses);
  }

  t8_cmesh_join_by_vertices_impl (cmesh, ntrees, eclasses, global_ids);
}

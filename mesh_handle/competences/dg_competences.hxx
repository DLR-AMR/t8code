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

/** \file dg_element_competences.hxx
 * TODO: IDEA for test: write rank on element data and then check if this is correct for the ghosts
 */

#pragma once

#include <t8.h>
#include <t8_types/t8_operators.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>
#include <optional>

namespace t8_mesh_handle
{
/// Dummy value for the local rank.
inline constexpr int LOCAL_RANK = -1;

/**
 * This mesh competence provides the competence to store a vector with the ranks of the according elements. 
 * This is always \a LOCAL_RANK for the local elements. The interesting parts are the ghost elements. 
 * The order of the vector is the same as in the mesh handle.
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct remote_ranks_mesh_competence: public t8_crtp_operator<TUnderlying, remote_ranks_mesh_competence>
{
 public:
  /**
   * Fill \a m_ranks with LOCAL_RANK for each local element and with the
   * corresponding remote MPI rank for each ghost element, i.e.:
   *   m_ranks[0 .. num_local_elements - 1]  = LOCAL_RANK
   *   m_ranks[num_local_elements .. num_local_elements + num_ghosts - 1]
   *                                         = remote rank of that ghost
   */
  void
  set_rank_vector () const
  {
    //TODO if no ghost is defined then return vec with only llocal ranks (num ghosts is 0 )
    const TUnderlying& mesh = this->underlying ();
    t8_forest_t forest = mesh.get_forest ();

    const t8_locidx_t num_local = mesh.get_num_local_elements ();
    const t8_locidx_t num_ghosts = mesh.get_num_ghosts ();

    m_ranks.resize (num_local + num_ghosts);

    /* --- local elements: rank = LOCAL_RANK --- */
    for (t8_locidx_t ilocal = 0; ilocal < num_local; ++ilocal) {
      m_ranks[ilocal] = LOCAL_RANK;
    }

    if (num_ghosts == 0) {
      return;
    }

    /* --- ghost elements: determine rank per ghost element --- */

    /* t8_forest_ghost_get_remotes returns a plain int array of length
     * num_remotes in ascending order, and fills num_remotes. */
    int num_remotes = 0;
    int* remotes = t8_forest_ghost_get_remotes (forest, &num_remotes);

    for (int iremote = 0; iremote < num_remotes; ++iremote) {
      const int remote = remotes[iremote];

      /* First ghost-element index (0-based within ghost elements only)
       * that belongs to this remote rank. */
      const t8_locidx_t first_elem = t8_forest_ghost_remote_first_elem (forest, remote);

      /* The number of ghost elements of this remote is the difference
       * between its first element index and the next remote's first
       * element index (or the total ghost count for the last remote). */
      t8_locidx_t num_elems_of_remote;
      if (iremote + 1 < num_remotes) {
        const int next_remote = remotes[iremote + 1];
        const t8_locidx_t next_first = t8_forest_ghost_remote_first_elem (forest, next_remote);
        num_elems_of_remote = next_first - first_elem;
      }
      else {
        num_elems_of_remote = num_ghosts - first_elem;
      }

      /* Write the remote rank for every ghost element of this remote.
       * The handle id of a ghost element is num_local + (ghost-only index). */
      for (t8_locidx_t ielem = 0; ielem < num_elems_of_remote; ++ielem) {
        m_ranks[num_local + first_elem + ielem] = remote;
      }
    }
  }

  int
  get_rank (t8_locidx_t element_handle_id)
  {
    T8_ASSERT (element_handle_id < m_ranks.size ());
    return m_ranks[element_handle_id];
  }

 protected:
  mutable std::vector<int> m_ranks;
};

struct face_side
{
  int m_element_id;         // local element id, or ghost layer index if remote
  int m_local_face_id;      // which face of that element points to this interface
  int m_orientation;        // face orientation code for coordinate permutation
  int m_rank = LOCAL_RANK;  // LOCAL_RANK if owned locally, else MPI rank of owner

  bool
  is_remote () const
  {
    return m_rank != LOCAL_RANK;
  }
};

enum class face_type {
  CONFORMAL,      // exactly 2 sides, same level, all local
  MORTAR,         // 1 large side + N small sides (N=2 in 2D, up to 4 in 3D hex)
  BOUNDARY,       // exactly 1 side, domain boundary
  MPI_CONFORMAL,  // exactly 2 sides, same level, at least one remote
  MPI_MORTAR,     // mortar where at least one side is remote
};

struct face
{
  face_type m_type;
  std::vector<face_side> m_sides;
  // Convention for MORTAR / MPI_MORTAR:
  //   m_sides[0]     = large side
  //   m_sides[1..N]  = small sides (in face-corner order of the large element)
  // Convention for CONFORMAL / MPI_CONFORMAL:
  //   m_sides[0]     = primary (lower global id or local owner)
  //   m_sides[1]     = secondary
  // Convention for BOUNDARY:
  //   m_sides[0]     = the single local side
  /** \note Currently unused but may be filled with t8_cmesh_set_attribute. */
  int m_boundary_tag = 0;  // meaningful only for BOUNDARY
};

/**
 * TODO
 * \tparam TUnderlying Use the \ref mesh with specified competences as template parameter.
 */
template <typename TUnderlying>
struct face_vector_mesh_competence: public t8_crtp_operator<TUnderlying, face_vector_mesh_competence>
{
 public:
  /**
   * Build \a m_faces so that every face touching at least one local element
   * appears exactly once.  Simultaneously populate \a m_element_face_vector
   * so that m_element_face_vector[handle_id] holds one index into m_faces
   * per face of that element.
   *
   * Uniqueness rules
   * ----------------
   * BOUNDARY      – only one local side exists, always added.
   * CONFORMAL /
   * MPI_CONFORMAL – added by the side whose element has the smaller
   *                 handle_id.  For a local–ghost pair the local element
   *                 always has the smaller id (local ids < ghost ids), so
   *                 the local side is always the one that inserts the face.
   * MORTAR /
   * MPI_MORTAR    – the large (coarser) side owns the face and inserts it.
   *                 The small sides are added to m_faces[*].m_sides by the
   *                 large-side owner; they do not trigger their own insertion.
   */
  void
  set_unique_face_vector () const
  {
    const TUnderlying& mesh = this->underlying ();
    if constexpr (requires (TUnderlying& meshvar) { meshvar.set_rank_vector (); })
      mesh.set_rank_vector ();
  }
  else
  {
    SC_ABORT ("Use remote_ranks_mesh_competence together with face_vector_mesh_competence.\n");
  }

  /* One entry per element (local + ghost), each entry is a list of
     * indices into m_faces. */
  m_element_face_vector.assign (mesh.get_num_local_elements () + mesh.get_num_ghosts (), {});

  for (const auto& elem : mesh) {
    const int handle_id = elem.get_element_handle_id ();

    for (int iface = 0; iface < elem.get_num_faces (); ++iface) {
      std::vector<int> dual_faces;
      auto neighs = elem.get_face_neighbors (iface, dual_faces);
      const int num_neighs = static_cast<int> (neighs.size ());

      if (num_neighs == 0) {
        /* ---- BOUNDARY ---- */
        face f;
        f.m_type = face_type::BOUNDARY;
        f.m_sides.push_back ({ handle_id, iface, /*orientation=*/0, LOCAL_RANK });

        const int face_idx = static_cast<int> (m_faces.size ());
        m_faces.push_back (std::move (f));
        m_element_face_vector[handle_id].push_back (face_idx);
      }
      else if (num_neighs == 1) {
        const int neigh_handle = neighs[0]->get_element_handle_id ();
        const int neigh_rank = mesh.get_rank (neigh_handle);
        const bool neigh_is_remote = (neigh_rank != LOCAL_RANK);

        if (elem.get_level () == neighs[0]->get_level ()) {
          /* ---- CONFORMAL or MPI_CONFORMAL ----
             * Insert only from the side with the smaller handle_id so that
             * each conformal face appears exactly once.  For a local–ghost
             * pair the local element always wins (local ids < ghost ids). */
          if (handle_id < neigh_handle) {
            face f;
            f.m_type = neigh_is_remote ? face_type::MPI_CONFORMAL : face_type::CONFORMAL;
            f.m_sides.push_back ({ handle_id, iface, /*orientation=*/0, LOCAL_RANK });
            f.m_sides.push_back ({ neigh_handle, dual_faces[0], /*orientation=*/0, neigh_rank });

            const int face_idx = static_cast<int> (m_faces.size ());
            m_faces.push_back (std::move (f));
            m_element_face_vector[handle_id].push_back (face_idx);
            m_element_face_vector[neigh_handle].push_back (face_idx);
          }
          else {
            /* The neighbour (local) will insert this face; just record
               * that we will participate so the index can be filled in when
               * the neighbour processes its face.  Because we iterate in
               * handle_id order the neighbour has already done so. */
            /* Find the face index already inserted by the neighbour. */
            for (int fi : m_element_face_vector[neigh_handle]) {
              const face& candidate = m_faces[fi];
              if (candidate.m_sides.size () >= 2
                  && ((candidate.m_sides[0].m_element_id == neigh_handle
                       && candidate.m_sides[1].m_element_id == handle_id)
                      || (candidate.m_sides[1].m_element_id == neigh_handle
                          && candidate.m_sides[0].m_element_id == handle_id))) {
                m_element_face_vector[handle_id].push_back (fi);
                break;
              }
            }
          }
        }
        else if (elem.get_level () > neighs[0]->get_level ()) {
          /* ---- Small side of a MORTAR or MPI_MORTAR ----
             * The large (coarser) neighbour owns and inserts this face.
             * We are on a small side; do not insert, but record the
             * face index once the large side has inserted it.
             * Because we may reach the small side before the large side
             * has been processed (it is a ghost), we defer: the large
             * side's insertion loop will also update our entry. */
          /* Nothing to insert here; the large side handles it. */
        }
        else {
          /* ---- Large side of a MORTAR or MPI_MORTAR (num_neighs==1) ----
             * Our level < neighbour level: we are the coarse side.
             * This branch is reached when a coarse local element borders a
             * single finer element (possible in 2D with triangles). */
          face f;
          f.m_type = neigh_is_remote ? face_type::MPI_MORTAR : face_type::MORTAR;
          f.m_sides.push_back ({ handle_id, iface, /*orientation=*/0, LOCAL_RANK });
          f.m_sides.push_back ({ neigh_handle, dual_faces[0], /*orientation=*/0, neigh_rank });

          const int face_idx = static_cast<int> (m_faces.size ());
          m_faces.push_back (std::move (f));
          m_element_face_vector[handle_id].push_back (face_idx);
          m_element_face_vector[neigh_handle].push_back (face_idx);
        }
      }
      else {
        /* ---- Large side of a MORTAR or MPI_MORTAR (num_neighs > 1) ----
           * We are the coarse element; neighs are all finer small sides. */
        bool any_remote = false;
        for (int in = 0; in < num_neighs; ++in) {
          if (neighs[in]->get_rank () != LOCAL_RANK) {
            any_remote = true;
            break;
          }
        }

        face f;
        f.m_type = any_remote ? face_type::MPI_MORTAR : face_type::MORTAR;
        /* m_sides[0] = large side (us) */
        f.m_sides.push_back ({ handle_id, iface, /*orientation=*/0, LOCAL_RANK });
        /* m_sides[1..N] = small sides in face-corner order */
        for (int in = 0; in < num_neighs; ++in) {
          const int nh = neighs[in]->get_element_handle_id ();
          const int nr = neighs[in]->get_rank ();
          f.m_sides.push_back ({ nh, dual_faces[in], /*orientation=*/0, nr });
        }

        const int face_idx = static_cast<int> (m_faces.size ());
        m_faces.push_back (std::move (f));

        /* Record this face index for the large side and all small sides. */
        m_element_face_vector[handle_id].push_back (face_idx);
        for (int in = 0; in < num_neighs; ++in) {
          m_element_face_vector[neighs[in]->get_element_handle_id ()].push_back (face_idx);
        }
      }
    }
  }
}

protected:
  mutable std::vector<face>
    m_faces;
/* For each element (local + ghost), the list of indices into m_faces
   * that this element participates in.  One entry per face of the element. */
mutable std::vector<std::vector<int>> m_element_face_vector;
};
}  // namespace t8_mesh_handle

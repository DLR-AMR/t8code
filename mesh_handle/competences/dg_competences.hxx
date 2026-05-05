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
 * TODO
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
 * \tparam TUnderlying Use the \ref element with specified competences as template parameter.
 */
template <typename TUnderlying>
struct face_vector_mesh_competence: public t8_crtp_operator<TUnderlying, face_vector_mesh_competence>
{
 public:
  /**
   * TODO
   */
  void
  set_unique_face_vector () const
  {
    for (const auto& elem : this->underlying ()) {
      for (int iface = 0; iface < elem.get_num_faces (); ++iface) {
        std::vector<int> dual_faces;
        auto neighs = elem.get_face_neighbors (iface, dual_faces);
        auto num_neighs = neighs.size ();
        if (num_neighs == 0) {
          // Face is a geometry boundary.
          m_faces.push_back ({ mortar_type::BOUNDARY, { elem.get_element_handle_id () } });
        }
        else if (num_neighs > 1) {
          // Face is a big mortar.
          std::vector<int> temp;
          temp.push_back (elem.get_element_handle_id ());
          for (int ineigh = 0; ineigh < num_neighs; ++ineigh) {
            temp.push_back (neighs[ineigh]->get_element_handle_id ());
          }
          m_faces.push_back ({ mortar_type::BIGMORTAR, temp });
        }
        else {  //num_neighs==1
          if (neighs[0]->get_level () == elem.get_level ()) {
            // For conformal faces, the face should only be added once to the unique faces vector.
            if (elem.get_element_handle_id () < neighs[0]->get_element_handle_id ()) {
              m_faces.push_back (
                { mortar_type::CONFORMAL, { elem.get_element_handle_id (), neighs[0]->get_element_handle_id () } });
            }
          }
          else {
            // Face is the small mortar.
            m_faces.push_back (
              { mortar_type::SMALLMORTAR, { elem.get_element_handle_id (), neighs[0]->get_element_handle_id () } });
          }
        }
      }
    }
  }

 protected:
  mutable std::vector<face> m_faces;
  // Inverse map: for each element, which face entries does it participate in?
  std::vector<std::vector<&face>> m_element_face_vector;
};
}  // namespace t8_mesh_handle

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

/** \file dg_competences.hxx
 * Mesh competences that are useful for discontinuous Galerkin methods. This includes the storage of the remote 
 * ranks of the ghost elements and the face connectivity.
 * TODO: IDEA for test: write rank on element data and then check if this is correct for the ghosts
 */

#pragma once

#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8.h>
#include <t8_types/t8_operators.hxx>
#include <t8_types/t8_vec.hxx>
#include <vector>
#include <algorithm>

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
   * Fill \a ranks with LOCAL_RANK for each local element and with the
   * corresponding remote MPI rank for each ghost element, i.e.:
   *   ranks[0 .. num_local_elements - 1]  = LOCAL_RANK
   *   ranks[num_local_elements .. num_local_elements + num_ghosts - 1]
   *                                         = remote rank of that ghost
   */
  void
  set_rank_vector () const
  {
    const TUnderlying& mesh = this->underlying ();
    const t8_locidx_t num_local = mesh.get_num_local_elements ();
    const t8_locidx_t num_ghosts = mesh.get_num_ghosts ();
    if (num_ghosts == 0) {
      ranks.assign (num_local, LOCAL_RANK);
      return;
    }

    /* --- local elements: rank = LOCAL_RANK --- */
    ranks.resize (num_local + num_ghosts);
    std::fill_n (ranks.begin (), num_local, LOCAL_RANK);

    /* --- ghost elements: determine rank per ghost element. --- */
    /* Get array with remote ranks in ascending order. */
    const t8_forest_t forest = mesh.get_forest ();
    int num_remotes = 0;
    int* remotes = t8_forest_ghost_get_remotes (forest, &num_remotes);

    for (int iremote = 0; iremote < num_remotes; ++iremote) {
      const int remote = remotes[iremote];
      /* First ghost-element index that belongs to this remote rank. */
      const t8_locidx_t first_elem = t8_forest_ghost_remote_first_elem (forest, remote);
      /* The number of ghost elements of this remote is the difference between first_elem and the next remote's first
       * element index. */
      t8_locidx_t num_elems_of_remote;
      if (iremote + 1 < num_remotes) {
        const int next_remote = remotes[iremote + 1];
        const t8_locidx_t next_first = t8_forest_ghost_remote_first_elem (forest, next_remote);
        num_elems_of_remote = next_first - first_elem;
      }
      else {
        num_elems_of_remote = num_ghosts - first_elem;
      }

      /* Write remote rank for every ghost element of this remote. */
      for (t8_locidx_t ielem = 0; ielem < num_elems_of_remote; ++ielem) {
        ranks[num_local + first_elem + ielem] = remote;
      }
    }
  }
  /**
   * Get the rank associated with a given element handle.
   * \param[in] element_handle_id The ID of the element handle.
   * \return The rank of the specified element if \ref set_rank_vector was called beforehand.
   */
  int
  get_rank (t8_locidx_t element_handle_id) const
  {
    T8_ASSERT (static_cast<size_t> (element_handle_id) < ranks.size ());
    return ranks[element_handle_id];
  }

 protected:
  mutable std::vector<int> ranks;  ///< The rank of the owner for each element.
};

/** Type of the face. */
enum class face_type {
  BOUNDARY,       ///< Exactly 1 side, domain boundary.
  CONFORMAL,      ///< Exactly 2 sides, same level, all local.
  MORTAR,         ///< 1 large side + N small sides, all local.
  MPI_CONFORMAL,  ///< Exactly 2 sides, same level, exactly one remote.
  MPI_MORTAR,     ///< Mortar where at least one side (large or small) is remote.
};

/** Class for the face side of an element. One \ref face can have multiple face sides of different elements. */
struct face_side
{
  int element_id;         ///< Mesh element handle id of the side's element.
  int local_face_id;      ///< The face of that element pointing to this interface.
  int orientation = 0;    ///< Face orientation code for coordinate permutation.
  int rank = LOCAL_RANK;  ///< LOCAL_RANK if owned locally, else MPI rank of owner.
};

/** Class for a face. A face can have multiple \ref face_side s if different elements faces share the same face.
 */
struct face
{
  face_type type;  ///< The face type of the face, see \ref face_type for the possible types and their meaning.
  /** The face sides of different elements adjacent to this face. 
   * Order conventions for the sides in the vector depend on the face type:
   * - BOUNDARY: sides[0]     = the single local side
   * - CONFORMAL / MPI_CONFORMAL: sides[0] = primary (smaller handle id); sides[1] = secondary.
   *                              For MPI_CONFORMAL the local side is always the primary side with the smaller handle id (local ids < ghost ids).
   * - MORTAR / MPI_MORTAR: sides[0] = large side; 
   *                        sides[1..N] = small sides (in face-corner order of the large element)
   */
  std::vector<face_side> sides;
};

/**
 * Mesh competence to build a unique vector of faces, where each face appears exactly once. 
 * This is useful for DG methods, where we need to loop over faces and access the neighboring elements across the face.
 * This mesh competence should be combined with the competence \ref remote_ranks_mesh_competence.
 * Each vector entry is of type \ref face and we have the following uniqueness rules for the face types specified in
 * \ref face_type :
 * - BOUNDARY: only one local side exists, always added.
 * - CONFORMAL / MPI_CONFORMAL: added by the side whose element has the smaller handle_id. 
 *      For a local–ghost pair the local element always has the smaller id (local ids < ghost ids),
 *      so the local side is always the one that inserts the face.
 * - MORTAR / MPI_MORTAR: the large (coarser) side owns the face and inserts it (also for ghosts). 
 *      The small sides are specified in sides.
 * Additionally, a vector is build that holds the face indices for each element.
 * 
 * \tparam TUnderlying Use the \ref mesh with specified competences as template parameter.
 */
template <typename TUnderlying>
struct face_vector_mesh_competence: public t8_crtp_operator<TUnderlying, face_vector_mesh_competence>
{
 public:
  /**
   * Build \a m_faces so that every face touching at least one local element appears exactly once. 
   * Simultaneously populate \a m_element_face_vector so that m_element_face_vector[handle_id] holds one index 
   * into m_faces per face of that element.
   */
  void
  set_unique_face_vector () const
  {
    if (!m_faces.empty () || !m_element_face_vector.empty ()) {
      m_faces.clear ();
      m_element_face_vector.clear ();
    }
    const TUnderlying& mesh = this->underlying ();
    T8_ASSERT (mesh.is_balanced ());
    // Ensure that rank vector is set.
    if constexpr (TUnderlying::has_remote_ranks_mesh_competence ()) {
      mesh.set_rank_vector ();
    }
    else {
      SC_ABORT ("Use remote_ranks_mesh_competence together with face_vector_mesh_competence.\n");
    }

    const t8_forest_t forest = mesh.get_forest ();
    SC_CHECK_ABORT (forest->incomplete_trees == 0, "This functionality is not specified for incomplete trees.\n");

    // Vector should store one entry per face of each element, so reserve space accordingly.
    m_element_face_vector.assign (mesh.get_num_local_elements () + mesh.get_num_ghosts (), {});

    // Check each face of each element and insert into vectors.
    for (const auto& elem : mesh) {
      const int handle_id = elem.get_element_handle_id ();
      const int num_faces = elem.get_num_faces ();
      if (m_element_face_vector[handle_id].empty ()) {
        m_element_face_vector[handle_id].assign (num_faces, -1);
      }
      for (int iface = 0; iface < num_faces; ++iface) {
        // Continue if the face entry is already filled.
        if (m_element_face_vector[handle_id][iface] != -1) {
          continue;
        }
        // Get neighbors and their dual face numbers to check the face_type.
        std::vector<int> dual_faces;
        auto neighs = elem.get_face_neighbors (iface, dual_faces);
        const int num_neighs = static_cast<int> (neighs.size ());

        if (num_neighs == 0) {
          /* --- BOUNDARY --- */
          face f { face_type::BOUNDARY, { { handle_id, iface, /*orientation=*/0, LOCAL_RANK } } };
          const int face_idx = static_cast<int> (m_faces.size ());
          m_faces.push_back (std::move (f));
          m_element_face_vector[handle_id][iface] = face_idx;
        }
        else if (num_neighs == 1) {
          const int neigh_id = neighs[0]->get_element_handle_id ();
          const int neigh_rank = mesh.get_rank (neigh_id);

          if (elem.get_level () == neighs[0]->get_level ()) {
            /* --- CONFORMAL or MPI_CONFORMAL --- */
            // Insert only from the side with the smaller handle_id so that each conformal face appears exactly once.
            if (handle_id < neigh_id) {
              face f;
              f.type = (neigh_rank != LOCAL_RANK) ? face_type::MPI_CONFORMAL : face_type::CONFORMAL;

              const int orientation = t8_forest_leaf_face_orientation (
                forest, elem.get_local_tree_id (), t8_forest_get_scheme (forest), elem.get_forest_element (), iface);
              f.sides.push_back ({ handle_id, iface, orientation, LOCAL_RANK });

              const int orientation_neigh = t8_forest_leaf_face_orientation (
                forest, neighs[0]->get_local_tree_id (), t8_forest_get_scheme (forest),
                neighs[0]->get_forest_element (), dual_faces[0]);
              f.sides.push_back ({ neigh_id, dual_faces[0], orientation_neigh, neigh_rank });

              const int face_idx = static_cast<int> (m_faces.size ());
              m_faces.push_back (std::move (f));
              m_element_face_vector[handle_id][iface] = face_idx;
              if (m_element_face_vector[neigh_id].empty ()) {
                m_element_face_vector[neigh_id].assign (neighs[0]->get_num_faces (), -1);
              }
              m_element_face_vector[neigh_id][dual_faces[0]] = face_idx;
            }
          }
          else if (elem.get_level () > neighs[0]->get_level ()) {
            /* --- Small side of a MORTAR or MPI_MORTAR --- */
            /* The large (coarser) neighbour owns and inserts this face.
             * We only need to do something here if the neighbor is remote.
             */
            if (neigh_rank == LOCAL_RANK) {
              continue;  // local large side will insert the face, so we can skip this.
            }

            // Check if the ghost is already inserted by another small mortar. If yes, only add the face.
            if (!m_element_face_vector[neigh_id].empty () && m_element_face_vector[neigh_id][dual_faces[0]] != -1) {
              const int orientation = t8_forest_leaf_face_orientation (
                forest, elem.get_local_tree_id (), t8_forest_get_scheme (forest), elem.get_forest_element (), iface);
              m_faces[m_element_face_vector[neigh_id][dual_faces[0]]].sides.push_back (
                { handle_id, iface, orientation, LOCAL_RANK });
              m_element_face_vector[handle_id][iface] = m_element_face_vector[neigh_id][dual_faces[0]];
              continue;
            }
            // Construct new face.
            const int orientation_neigh
              = t8_forest_leaf_face_orientation (forest, neighs[0]->get_local_tree_id (), t8_forest_get_scheme (forest),
                                                 neighs[0]->get_forest_element (), dual_faces[0]);
            const int orientation = t8_forest_leaf_face_orientation (
              forest, elem.get_local_tree_id (), t8_forest_get_scheme (forest), elem.get_forest_element (), iface);
            face f { face_type::MPI_MORTAR,
                     { { neigh_id, dual_faces[0], orientation_neigh, neigh_rank },  // Large mortar first.
                       { handle_id, iface, orientation, LOCAL_RANK } } };

            const int face_idx = static_cast<int> (m_faces.size ());
            m_faces.push_back (std::move (f));
            m_element_face_vector[handle_id][iface] = face_idx;
            if (m_element_face_vector[neigh_id].empty ()) {
              m_element_face_vector[neigh_id].assign (neighs[0]->get_num_faces (), -1);
            }
            m_element_face_vector[neigh_id][dual_faces[0]] = face_idx;
          }
          else {
            SC_ABORT ("Not possible without incomplete trees.\n");
          }
        }
        else {
          /* --- Large side of a MORTAR or MPI_MORTAR (num_neighs > 1) ---  */
          // We are the coarse element; neighs are all small mortars.
          bool any_remote = std::any_of (neighs.begin (), neighs.end (), [&] (const auto* neigh) {
            return mesh.get_rank (neigh->get_element_handle_id ()) != LOCAL_RANK;
          });
          // Face index for the face we are about to insert.
          const int face_idx = static_cast<int> (m_faces.size ());
          m_element_face_vector[handle_id][iface] = face_idx;

          face f;
          f.type = any_remote ? face_type::MPI_MORTAR : face_type::MORTAR;
          const int orientation = t8_forest_leaf_face_orientation (
            forest, elem.get_local_tree_id (), t8_forest_get_scheme (forest), elem.get_forest_element (), iface);
          f.sides.push_back ({ handle_id, iface, orientation, LOCAL_RANK });
          for (int ineigh = 0; ineigh < num_neighs; ++ineigh) {
            const int neigh_id = neighs[ineigh]->get_element_handle_id ();
            const int neigh_rank = mesh.get_rank (neigh_id);
            const int orientation_neigh = t8_forest_leaf_face_orientation (
              forest, neighs[ineigh]->get_local_tree_id (), t8_forest_get_scheme (forest),
              neighs[ineigh]->get_forest_element (), dual_faces[ineigh]);
            f.sides.push_back ({ neigh_id, dual_faces[ineigh], orientation_neigh, neigh_rank });
            // Record face index for the small mortar.
            if (m_element_face_vector[neigh_id].empty ()) {
              m_element_face_vector[neigh_id].assign (neighs[ineigh]->get_num_faces (), -1);
            }
            m_element_face_vector[neigh_id][dual_faces[ineigh]] = face_idx;
          }
          // Insert face f.
          m_faces.push_back (std::move (f));
        }
      }
    }
  }

  /** Get the vector of unique faces. \ref set_unique_face_vector must be called beforehand to populate the vector.
   * \return Const reference to the vector of faces.
   */
  const std::vector<face>&
  get_unique_face_vector () const
  {
    T8_ASSERTF (!m_faces.empty (), "m_faces has not been set. Call set_unique_face_vector() first.");
    return m_faces;
  }

  /** Get the vector with indices of faces for each mesh handle element.
 * \ref set_unique_face_vector must be called beforehand to populate the vector.
 * \return Const reference to the element-face vector.
 */
  const std::vector<std::vector<int>>&
  get_element_face_vector () const
  {
    T8_ASSERTF (!m_element_face_vector.empty (),
                " m_element_face_vector has not been set. Call set_unique_face_vector() first.");
    return m_element_face_vector;
  }

 protected:
  /** Vector with one entry per unique face. */
  mutable std::vector<face> m_faces;
  /** For each element (local + ghost), the list of indices into m_faces
   * that this element participates in. One entry per face of the element. */
  mutable std::vector<std::vector<int>> m_element_face_vector;
};
}  // namespace t8_mesh_handle

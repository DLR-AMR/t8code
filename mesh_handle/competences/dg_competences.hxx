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
// Sentinel value for "this side is local"
inline constexpr int LOCAL_RANK = -1;

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
  int m_boundary_tag = 0;  // meaningful only for BOUNDARY // \note
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
    for (const auto &elem : this->underlying ()) {
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
  // Store as (type, index_into_vector, side_role) for each local face
  std::vector<std::array<std::optional<face_ref>, MAX_FACES_PER_ELEM>> m_element_face_map;
};
}  // namespace t8_mesh_handle

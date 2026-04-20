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

struct mortar_type
{
  static const int CONFORMAL = 0;
  static const int SMALLMORTAR = 1;
  static const int BIGMORTAR = 2;
  static const int BOUNDARY = -1;
} struct face
{
  mortar_type m_mortartype;
  std::vector<int> m_element_ids;
  std::optional<int> ranks;
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
        auto neighs = elem.get_face_neighbors (iface);
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
  mutable std::vector<std::vector<const face &>> m_elements_to_faces;
};
}  // namespace t8_mesh_handle

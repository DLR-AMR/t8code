/*
This file is part of t8code.
t8code is a C library to manage a collection (a forest) of multiple
connected adaptive space-trees of general element classes in parallel.

Copyright (C) 2025 the developers

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

/** \file t8_unstructured_ghost.hxx
 * TODO
 */

#ifndef T8_UNSTRUCTURED_GHOST_HXX
#define T8_UNSTRUCTURED_GHOST_HXX

#include <t8.h>
#include <t8_unstructured_mesh/t8_unstructured_element.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <vector>

template <template <typename> class... TCompetence>
class t8_unstructured_ghost_element: public t8_unstructured_mesh_element<TCompetence...> {
  using Base = t8_unstructured_mesh_element<TCompetence...>;

  t8_unstructured_ghost_element (t8_unstructured_mesh<TCompetence...>* unstructured_mesh, t8_locidx_t tree_id,
                                 t8_locidx_t lghost_tree_id, t8_locidx_t element_id)
    : Base (unstructured_mesh, tree_id, element_id), m_lghost_tree_id (lghost_tree_id)
  {
  }

  std::vector<t8_locidx_t>
  get_face_neighbors (int face, int* num_neighbors, int* dual_faces[]) = delete;

 private:
  const t8_element_t*
  get_element () const
  {
    return t8_forest_ghost_get_leaf_element (m_unstructured_mesh->m_forest, m_lghost_tree_id, m_element_id);
  }

  t8_locidx_t m_lghost_tree_id;
};

#endif

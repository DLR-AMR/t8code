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

/** \file ghost.hxx
 * TODO
 */

#ifndef T8_GHOST_HXX
#define T8_GHOST_HXX

#include <t8.h>
#include <t8_mesh_handle/element.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <vector>

namespace t8_mesh_handle
{
/* Forward declaration of the mesh class of the handle.
 */
template <template <typename> class... TCompetence>
class mesh;

template <template <typename> class... TCompetence>
class ghost_element: public element<TCompetence...> {
  using Base = element<TCompetence...>;
  using mesh_class = mesh<TCompetence...>;
  friend mesh_class;

 protected:
  ghost_element (mesh<TCompetence...>* mesh, t8_locidx_t tree_id, t8_locidx_t lghost_tree_id, t8_locidx_t element_id)
    : Base (mesh, tree_id, element_id), m_lghost_tree_id (lghost_tree_id)
  {
  }

 public:
  /**
   * TODO
   */
  static constexpr bool
  is_ghost_element ()
  {
    return true;
  }

  /** It is not possible for ghosts to compute their neighbors. */
  std::vector<t8_locidx_t>
  get_face_neighbors (int face, int* num_neighbors, int* dual_faces[]) = delete;

 private:
  const t8_element_t*
  get_element () const
  {
    return t8_forest_ghost_get_leaf_element (this->m_mesh->m_forest, m_lghost_tree_id, this->m_element_id);
  }

  t8_locidx_t m_lghost_tree_id;
};

}  // namespace t8_mesh_handle
#endif /* !T8_GHOST_HXX */

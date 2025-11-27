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

/** \file ghost_element.hxx
 * Implementation of the ghost element class of the \ref t8_mesh_handle::mesh handle.
 */

#ifndef T8_GHOST_ELEMENT_HXX
#define T8_GHOST_ELEMENT_HXX

#include <t8.h>
#include "abstract_element.hxx"
#include <t8_forest/t8_forest_ghost.h>

namespace t8_mesh_handle
{
/* Forward declaration of the \ref mesh class of the handle.
 */
template <template <typename> class... TCompetence>
class mesh;

/** 
 * Class for the ghost elements of the \ref mesh handle. 
 * This class is a child class of the abstract element class and implements the functionality specific to ghost elements.
 * See \ref abstract_element for more information, especially on the template parameter.
 *
 * \tparam TCompetence The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class ghost_element: public abstract_element<TCompetence...> {
  using Base = abstract_element<TCompetence...>;
  using mesh_class = mesh<TCompetence...>;
  friend mesh_class;

  /**
   * Constructor for a ghost element. This constructor can only be called by the friend mesh class as
   * the user should use the mesh class to construct ghost elements.
   * \param [in] mesh             Pointer to the mesh the ghost element should belong to.
   * \param [in] tree_id          The tree id of the element (normally the number of local trees + local ghost tree id).
   * \param [in] lghost_tree_id   The ghost tree id of the element in the forest defining the mesh.
   * \param [in] element_id       The element id of the element in the forest defining the mesh.
   */
  ghost_element (mesh_class* mesh, t8_locidx_t tree_id, t8_locidx_t lghost_tree_id, t8_locidx_t element_id)
    : Base (mesh, tree_id, element_id), m_lghost_tree_id (lghost_tree_id)
  {
  }

 public:
  /**
   * Implementation of the function to check if the element is a ghost element.
   * \return Always true for the ghost element class.
   */
  constexpr bool
  is_ghost_element () const override
  {
    return true;
  }

 private:
  /**
   * Implementation of the getter for the leaf element of the ghost element.
   * \return The leaf element.
   */
  const t8_element_t*
  get_element () const override
  {
    return t8_forest_ghost_get_leaf_element (this->m_mesh->m_forest, m_lghost_tree_id, this->m_element_id);
  }

  t8_locidx_t m_lghost_tree_id;
};

}  // namespace t8_mesh_handle
#endif /* !T8_GHOST_ELEMENT_HXX */

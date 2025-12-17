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

/** \file TODO
 */

#ifndef T8_MESH_HXX
#define T8_MESH_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include "mesh.hxx"
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <type_traits>
#include <functional>

namespace t8_mesh_handle
{

template <typename mesh_class>
/** Singleton for callback. Problem is the mesh class template!*/
class CallbackRegistry {
 public:
  using CoarsenCallBackType = std::function<bool (mesh_class mesh, std::vector<mesh_element>& family)>;

  static void
  register_coarsen_callback (CoarsenCallBackType Coarsencallback)
  {
    get_map ()[name] = std::move (cb);
  }

  static CallbackType
  get_coarsen_callback ()
  {
    auto& map = get_map ();
    auto it = map.find (name);
    if (it != map.end ()) {
      return it->second;
    }
    return nullptr;
  }

 private:
  static std::unordered_map<std::string, CallbackType>&
  get_map ()
  {
    static std::unordered_map<std::string, CallbackType> map;
    return map;
  }
};

template <typename mesh_class>
using coarsen_mesh_element_family = bool (*) (mesh_class mesh, std::vector<mesh_class::mesh_element_class>& family);

template <typename mesh_class>
using refine_mesh_element = bool (*) (mesh_class mesh, mesh_class::mesh_element_class& element);

template <typename mesh_class, typename coarsen_mesh_element_family, typename refine_mesh_element>
int
mesh_adapt_callback (mesh_class mesh, t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                     [[maybe_unused]] t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                     [[maybe_unused]] const t8_scheme* scheme, const int is_family,
                     [[maybe_unused]] const int num_elements, t8_element_t* elements[])
{

  if (refine_mesh_element (mesh, mesh.get_mesh_element (lelement_id))) {
    /* Refine this element. */
    return 1;
  }
  else if (is_family && coarsen_mesh_element_family (mesh, )) {
    /* Coarsen this family. Note that we check for is_family before, since returning < 0
     * if we do not have a family as input is illegal. */
    return -1;
  }
  /* Do not change this element. */
  return 0;
}

template <typename mesh_class>
void
adapt_mesh (mesh_class& mesh_handle,
            refine_mesh_element<mesh_class, typename mesh_class::mesh_element_class> refinement_callback,
            coarsen_mesh_element_family<mesh_class, typename mesh_class::mesh_element_class> coarsen_callback,
            bool recursive)
{
  auto forest_from = mesh_handle.get_forest ();

  t8_forest_t forest;
  t8_forest_init (&forest);
  t8_forest_set_adapt (forest, forest_from, mesh_adapt_callback<mesh_class, coarsen_callback, refinement_callback>,
                       recursive);
  t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  t8_forest_set_user_data (forest, mesh_handle.get_user_data ());

  t8_forest_commit (forest);
  mesh_handle.set_forest (forest);
}

}  // namespace t8_mesh_handle
#endif /* !T8_MESH_HXX */

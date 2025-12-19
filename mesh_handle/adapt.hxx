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

/** \file adapt.hxx
 */

#ifndef T8_ADAPT_HXX
#define T8_ADAPT_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include "mesh.hxx"
#include <t8_forest/t8_forest_ghost.h>
#include <vector>
#include <type_traits>
#include <functional>

namespace t8_mesh_handle
{

/** TODO*/
template <typename TMeshClass>
using coarsen_mesh_element_family
  = std::function<bool (TMeshClass, std::vector<typename TMeshClass::mesh_element_class>&)>;

template <typename TMeshClass>
using refine_mesh_element = std::function<bool (TMeshClass, typename TMeshClass::mesh_element_class&)>;

namespace detail
{

struct MeshAdaptContextBase
{
  virtual ~MeshAdaptContextBase () = default;

  virtual int
  adapt_callback (t8_locidx_t lelement_id, int is_family, int num_elements, t8_element_t* elements[])
    = 0;
};

template <typename TMesh>
struct MeshAdaptContext final: MeshAdaptContextBase
{
  using Element = typename TMesh::mesh_element_class;

  MeshAdaptContext (TMesh& m, refine_mesh_element<TMesh> r, coarsen_mesh_element_family<TMesh> c)
    : mesh (m), refine (std::move (r)), coarsen (std::move (c))
  {
  }

  int
  adapt_callback (t8_locidx_t lelement_id, int is_family, int num_elements, t8_element_t* elements[]) override
  {
    // refine
    if (refine) {
      Element elem = mesh.get_mesh_element (lelement_id);
      if (refine (mesh, elem)) {
        return 1;
      }
    }

    // coarsen
    if (is_family && coarsen) {
      std::vector<Element> element_family;
      for (int i = 0; i < num_elements; i++) {
        element_family.push_back (mesh.get_mesh_element (lelement_id + i));
      }
      if (coarsen (mesh, element_family)) {
        return -1;
      }
    }

    return 0;
  }

 private:
  TMesh& mesh;
  refine_mesh_element<TMesh> refine;
  coarsen_mesh_element_family<TMesh> coarsen;
};

/** Singleton for callback. Problem is the mesh class template!*/
class AdaptRegistry {
 public:
  static void
  register_context (t8_forest_t forest, MeshAdaptContextBase* context)
  {
    get_map ()[forest] = context;
  }

  static void
  unregister_context (t8_forest_t forest)
  {
    get_map ().erase (forest);
  }

  static MeshAdaptContextBase*
  get (t8_forest_t forest)
  {
    auto& map = get_map ();
    auto it = map.find (forest);
    return it != map.end () ? it->second : nullptr;
  }

 private:
  // Use getter instead of private member variable to ensure single initialization
  static std::unordered_map<t8_forest_t, MeshAdaptContextBase*>&
  get_map ()
  {
    static std::unordered_map<t8_forest_t, MeshAdaptContextBase*> map;
    return map;
  }
};

int
mesh_adapt_callback ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                     [[maybe_unused]] t8_eclass_t tree_class, t8_locidx_t lelement_id,
                     [[maybe_unused]] const t8_scheme* scheme, const int is_family, const int num_elements,
                     t8_element_t* elements[])
{
  auto* context = AdaptRegistry::get (forest_from);
  if (!context) {
    t8_global_infof (
      "Something went wrong while registering the adaption callbacks. Please check your implementation.");
    return 0;  // No adaption as default.
  }
  t8_locidx_t local_flat_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  // TODO: adapt names in other callbacks to not get confused between flat and this element_in_tree_index.
  return context->adapt_callback (local_flat_index, is_family, num_elements, elements);
}
}  // namespace detail

/** TODO Adapt a mesh handle according to the provided callbacks.
 * By default, the forest takes ownership of the source \b set_from such that it
 * will be destroyed on calling \ref t8_forest_commit. To keep ownership of \b
 * set_from, call \ref t8_forest_ref before passing it into this function.
 * This means that it is ILLEGAL to continue using \b set_from or dereferencing it
 * UNLESS it is referenced directly before passing it into this function.
 * \param [in,out] forest   The forest
 * \param [in] set_from     The source forest from which \b forest will be adapted.
 *                          We take ownership. This can be prevented by
 *                          referencing \b set_from.
 *                          If NULL, a previously (or later) set forest will
 *                          be taken (\ref t8_forest_set_partition, \ref t8_forest_set_balance).
 * \param [in] adapt_fn     The adapt function used on committing.
 * \param [in] recursive    A flag specifying whether adaptation is to be done recursively
 *                          or not. If the value is zero, adaptation is not recursive
 *                          and it is recursive otherwise.
 */
template <typename TMesh>
void
adapt_mesh (TMesh& mesh_handle, refine_mesh_element<TMesh> refine_callback,
            coarsen_mesh_element_family<TMesh> coarsen_callback, bool recursive)
{
  auto forest_from = mesh_handle.get_forest ();

  t8_forest_t forest;
  t8_forest_init (&forest);

  auto* context
    = new detail::MeshAdaptContext<TMesh> (mesh_handle, std::move (refine_callback), std::move (coarsen_callback));

  detail::AdaptRegistry::register_context (forest_from, context);

  t8_forest_set_adapt (forest, forest_from, detail::mesh_adapt_callback, recursive);
  t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
  t8_forest_set_user_data (forest, t8_forest_get_user_data (forest_from));

  t8_forest_commit (forest);
  mesh_handle.set_forest (forest);

  detail::AdaptRegistry::unregister_context (forest_from);
  delete context;
}

}  // namespace t8_mesh_handle
#endif /* !T8_ADAPT_HXX */

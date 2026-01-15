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
 * This file provides the functionality to adapt a \ref t8_mesh_handle::mesh
 * according to user defined callbacks.
 */

#ifndef T8_ADAPT_HXX
#define T8_ADAPT_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include "mesh.hxx"
#include <vector>
#include <functional>

namespace t8_mesh_handle
{

/** Namespace detail to hide implementation details from the user. */
namespace detail
{

/** Virtual base class for mesh adaptation contexts.
 * We need this base class and not only \ref MeshAdaptContext for the \ref AdaptRegistry. 
 * AdaptRegistry should not be templated because we need to access registered contexts in \ref mesh_adapt_callback_wrapper,
 * where we do not know the type of the mesh. Therefore, we work with a map of forests to instances of this base class to remain template free.
 */
struct MeshAdaptContextBase
{
  /** Virtual destructor for safe polymorphic deletion.
   */
  virtual ~MeshAdaptContextBase () = default;

  /** Pure virtual callback for mesh adaptation.
   * \param[in] lelement_handle_id Local element ID in the mesh handle.
   * \param [in] is_family         If 1, the entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements      The number of entries in \a elements.
   * \param [in] elements          Pointers to a family or, if \a is_family is zero, pointer to one element.
   *  \return 1 if the first entry in \a elements should be refined,
   *         -1 if the family \a elements shall be coarsened,
   *          0 else.
   */
  virtual int
  adapt_callback (const t8_locidx_t lelement_handle_id, const int is_family, const int num_elements,
                  t8_element_t* elements[])
    = 0;
};

/** Mesh adaptation context holding the mesh handle and the user defined callbacks.
 * Class inherits from \ref MeshAdaptContextBase and implements the virtual adapt callback using the mesh and the callbacks.
 * \tparam TMesh The mesh handle class.
 */
template <typename TMesh>
struct MeshAdaptContext final: MeshAdaptContextBase
{
  using Element = typename TMesh::element_class; /**< Type alias for the element class. */

  /** Constructor of the context with the mesh handle and the user defined callbacks.
   * \param [in] mesh_handle      The mesh handle to adapt.
   * \param [in] refine_callback  The refinement callback.
   * \param [in] coarsen_callback The coarsening callback.
   */
  MeshAdaptContext (TMesh& mesh_handle, refine_element<TMesh> refine_callback,
                    coarsen_element_family<TMesh> coarsen_callback)
    : m_mesh_handle (mesh_handle), m_refine_callback (std::move (refine_callback)),
      m_coarsen_callback (std::move (coarsen_callback))
  {
  }

  /** Callback for mesh adaptation using user defined callbacks.
   * \param [in] lelement_handle_id Local flat element ID in the mesh handle.
   * \param [in] is_family          If 1, the entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements       The number of entries in \a elements.
   * \param [in] elements           Pointers to a family or, if \a is_family is zero, pointer to one element.
   * \return 1 if the first entry in \a elements should be refined,
   *        -1 if the family \a elements shall be coarsened,
   *         0 else.
   */
  int
  adapt_callback (const t8_locidx_t lelement_handle_id, const int is_family, const int num_elements,
                  t8_element_t* elements[]) override
  {
    // Check if refine callback is set and call it using the correct mesh handle function arguments.
    if (m_refine_callback) {
      Element elem = m_mesh_handle[lelement_handle_id];
      if (m_refine_callback (m_mesh_handle, elem)) {
        return 1;
      }
    }

    // Check if a family is provided for adaption, if coarsen callback is set and call it using the correct mesh handle function arguments.
    if (is_family && m_coarsen_callback) {
      std::vector<Element> element_family;
      for (int i = 0; i < num_elements; i++) {
        element_family.push_back (m_mesh_handle[lelement_handle_id + i]);
      }
      if (m_coarsen_callback (m_mesh_handle, element_family)) {
        return -1;
      }
    }
    // Do nothing if the callbacks no refinement or coarsening should be done according to the callbacks.
    return 0;
  }

 private:
  TMesh& m_mesh_handle;                             /**< The mesh handle to adapt. */
  refine_element<TMesh> m_refine_callback;          /**< The refinement callback. */
  coarsen_element_family<TMesh> m_coarsen_callback; /**< The coarsening callback. */
};

/** Registry pattern is used to register contexts, which provides access to the callbacks and the mesh handle.
 * This globally accessible static class is required to get the handle and the callbacks in the forest callback, 
 * as the predefined header permits to give these as function arguments. 
 */
class AdaptRegistry {
 public:
  /** Static function to register \a context using \a forest as identifier. 
   * This makes the context publicly available using the Registry.
   * \param [in] forest The forest identifier.
   * \param [in] context The context to register.
   */
  static void
  register_context (t8_forest_t forest, MeshAdaptContextBase& context)
  {
    auto& map = get_map ();
    auto [it, inserted] = map.emplace (forest, context);
    if (!inserted) {
      throw std::logic_error ("Context already registered");
    }
  }

  /** Static function to unregister a context using \a forest as identifier. 
   * \param [in] forest The forest identifier.
   */
  static void
  unregister_context (t8_forest_t forest)
  {
    auto& map = get_map ();
    const auto erased = map.erase (forest);
    T8_ASSERT (erased == 1);
  }

  /** Getter for a context using \a forest as identifier. 
   * \param [in] forest The forest identifier.
   * \return Pointer to the context registered with the id \a forest if found, nullptr otherwise.
   */
  static MeshAdaptContextBase*
  get (t8_forest_t forest)
  {
    auto& map = get_map ();
    auto it = map.find (forest);
    return it != map.end () ? &it->second.get () : nullptr;
  }

 private:
  /** Get the static map associating t8_forest_t with MeshAdaptContextBase references.
   * We use a getter instead of private member variable to ensure single initialization
   * \return Reference to the static unordered map of t8_forest_t to MeshAdaptContextBase references.
   */
  static std::unordered_map<t8_forest_t, std::reference_wrapper<MeshAdaptContextBase>>&
  get_map ()
  {
    static std::unordered_map<t8_forest_t, std::reference_wrapper<MeshAdaptContextBase>> map;
    return map;
  }
};

/** Wrapper around the mesh handle adapt functionality to be able to pass the callbacks to the classic adapt routine of a forest. 
 * The function header fits the definition of \ref t8_forest_adapt_t.
 * \param [in] forest       Unused; forest to which the new elements belong.
 * \param [in] forest_from  Forest that is adapted.
 * \param [in] which_tree   Local tree containing \a elements.
 * \param [in] tree_class   Unused; eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in the tree of the first element in elements.
 * \param [in] scheme       Unused; scheme of the forest.
 * \param [in] is_family    If 1, the entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements.
 * \param [in] elements     Pointers to a family or, if \a is_family is zero, pointer to one element.
 * \return 1 if the first entry in \a elements should be refined,
 *        -1 if the family \a elements shall be coarsened,
 *         0 else.
 */
int
mesh_adapt_callback_wrapper ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             [[maybe_unused]] t8_eclass_t tree_class, t8_locidx_t lelement_id,
                             [[maybe_unused]] const t8_scheme* scheme, const int is_family, const int num_elements,
                             t8_element_t* elements[])
{
  // Get static adapt context from the registry.
  // Via this, we can access the mesh handle and the user defined callbacks that are using mesh handle functionality.
  auto* context = AdaptRegistry::get (forest_from);
  if (!context) {
    t8_global_infof (
      "Something went wrong while registering the adaption callbacks. Please check your implementation.");
    return 0;  // No adaption as default.
  }
  // Convert to index used in the mesh handle.
  t8_locidx_t mesh_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  // Call the actual adapt callback stored in the context.
  return context->adapt_callback (mesh_index, is_family, num_elements, elements);
}

}  // namespace detail

/** Adapt a mesh handle according to the provided callbacks.
 * The forest used to define the mesh handle is replaced in this function. 
 * If you want to additionally keep the current forest, call \ref t8_forest_ref before this function.
 * \param [in,out] mesh_handle   The mesh handle to adapt.
 * \param [in] refine_callback   Callback to decide if a single element should be refined.
 * \param [in] coarsen_callback  Callback to decide if a family of elements should be coarsened.
 * \param [in] recursive         Specifying whether adaptation is to be done recursively or not. 
 * \tparam TMesh                 The mesh handle class.
 */
// template <typename TMesh>
// void
// adapt_mesh (TMesh& mesh_handle, refine_element<TMesh> refine_callback, coarsen_element_family<TMesh> coarsen_callback,
//             bool recursive)
// {
//   auto forest_from = mesh_handle.get_forest ();
//   // Initialize forest for the adapted mesh.
//   t8_forest_t forest;
//   t8_forest_init (&forest);

//   // Create and register adaptation context holding the mesh handle and the user defined callbacks.
//   auto context
//     = detail::MeshAdaptContext<TMesh> (mesh_handle, std::move (refine_callback), std::move (coarsen_callback));
//   detail::AdaptRegistry::register_context (forest_from, context);

//   // Set up the forest for adaptation using the wrapper callback.
//   t8_forest_set_adapt (forest, forest_from, detail::mesh_adapt_callback_wrapper, recursive);
//   t8_forest_set_ghost (forest, 1, T8_GHOST_FACES);
//   t8_forest_set_user_data (forest, t8_forest_get_user_data (forest_from));
//   t8_forest_commit (forest);

//   // Replace the forest in the mesh handle with the adapted forest.
//   mesh_handle.set_forest (forest);
//   // Clean up.
//   detail::AdaptRegistry::unregister_context (forest_from);
// }

}  // namespace t8_mesh_handle
#endif /* !T8_ADAPT_HXX */

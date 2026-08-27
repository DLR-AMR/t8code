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

/** \file adapt.hxx
 * This file provides helper functionality to adapt a \ref t8_mesh_handle::mesh
 * according to a user defined callback.
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <mesh_handle/mesh.hxx>
#include <memory>
#include <span>

namespace t8_mesh_handle
{

/** Namespace detail to hide implementation details from the user. */
namespace detail
{

/** Virtual base class for mesh adaptation contexts.
 * We need this base class and not only \ref mesh_adapt_context for the \ref adapt_registry. 
 * adapt_registry should not be templated because we need to access registered contexts in \ref mesh_adapt_callback_wrapper,
 * where we do not know the type of the mesh. Therefore, we work with a map of forests to instances of this base class to remain template free.
 */
struct mesh_adapt_context_base
{
  /** Virtual destructor for safe polymorphic deletion.
   */
  virtual ~mesh_adapt_context_base () = default;

  /** Pure virtual callback for mesh adaptation.
   * \param [in] lelement_handle_id Local flat element ID in the mesh handle of the first element selected for adaptation.
   * \param [in] num_elements       The number of elements that should be considered. 
   *              If >1, the elements are expected to form a family that can be coarsened.
   * \return 1 if the element with index \a lelement_handle_id should be refined,
   *        -1 if the family shall be coarsened,
   *         0 else.
   */
  virtual int
  adapt_mesh (const t8_locidx_t lelement_handle_id, const int num_elements)
    = 0;
};

/** Templated mesh adaptation context holding the mesh handle and the user defined callback.
 * Struct inherits from \ref mesh_adapt_context_base and implements the virtual adapt callback using the mesh and the callback.
 * \tparam TMeshClass The mesh handle class.
 */
template <typename TMeshClass>
struct mesh_adapt_context final: mesh_adapt_context_base
{
  /** Constructor of the context with the mesh handle and the user defined callback.
   * \param [in] mesh_handle      The mesh handle to adapt.
   * \param [in] adapt_callback   The adapt callback.
   */
  mesh_adapt_context (TMeshClass& mesh_handle, typename TMeshClass::adapt_callback_type adapt_callback)
    : m_mesh_handle (mesh_handle), m_adapt_callback (std::move (adapt_callback))
  {
  }

  /** Callback for mesh adaptation using the user defined adapt callback.
   * \param [in] lelement_handle_id Local flat element ID in the mesh handle of the first element.
   * \param [in] num_elements       The number of elements that should be considered. 
   *              If >1, the elements are expected to form a family that can be coarsened.
   * \return 1 if the element with index \a lelement_handle_id should be refined,
   *        -1 if the family shall be coarsened,
   *         0 else.
   */
  int
  adapt_mesh (const t8_locidx_t lelement_handle_id, const int num_elements) override
  {
    // Check if adapt callback is set and call it using the correct mesh handle function arguments.
    T8_ASSERTF (m_adapt_callback, "No adapt callback set.");
    std::span<const typename TMeshClass::element_class> element_view (&m_mesh_handle[lelement_handle_id], num_elements);
    return m_adapt_callback (m_mesh_handle, element_view);
  }

 private:
  TMeshClass& m_mesh_handle;                                 /**< The mesh handle to adapt. */
  typename TMeshClass::adapt_callback_type m_adapt_callback; /**< The adapt callback. */
};

/** Registry pattern is used to register contexts, which provides access to the adapt callback and the mesh handle.
 * This globally accessible static class is required to get the handle and the callback in the forest callback, 
 * as the predefined header permits to give these as function arguments. 
 */
class adapt_registry {
 public:
  /** Static function to register \a context using \a forest as identifier. 
   * This makes the context publicly available using the registry.
   * \param [in] forest  The forest identifier. In our case, this is the forest from which we adapt the mesh because 
   *                     we do not change this forest during adaptation such that it is valid unique identifier.
   * \param [in] context The context to register. Use unique pointer to ensure proper memory management and ownership.
   */
  static void
  register_context (t8_forest_t forest, std::unique_ptr<mesh_adapt_context_base> context)
  {
    auto& map = get_map ();
    auto [it, inserted] = map.emplace (forest, std::move (context));
    if (!inserted) {
      t8_global_errorf ("ERROR: Context already registered!");
    }
  }

  /** Static function to unregister a context using \a forest as identifier. 
   * \param [in] forest The forest identifier. In our case, this is the forest from which we adapt.
   */
  static void
  unregister_context (t8_forest_t forest)
  {
    auto& map = get_map ();
    [[maybe_unused]] const auto erased = map.erase (forest);
    T8_ASSERT (erased == 1);
  }

  /** Getter for a context using \a forest as identifier. 
   * \param [in] forest The forest identifier. In our case, this is the forest from which we adapt.
   * \return Pointer to the context registered with the id \a forest if found, nullptr otherwise.
   */
  static mesh_adapt_context_base*
  get (t8_forest_t forest)
  {
    const auto& map = get_map ();
    const auto it = map.find (forest);
    return it != map.end () ? it->second.get () : nullptr;
  }

 private:
  /** Get the static map associating t8_forest_t with mesh_adapt_context_base references.
   * We use a getter instead of private member variable to ensure single initialization.
   * \return Reference to the static unordered map of t8_forest_t to mesh_adapt_context_base references.
   */
  static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_adapt_context_base>>&
  get_map ()
  {
    static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_adapt_context_base>> map;
    return map;
  }
};

/** Wrapper around the mesh handle adapt functionality to be able to pass the callback to the classic adapt routine of a forest. 
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
 * \note We do not support removing elements for the mesh handle.
 */
int
mesh_adapt_callback_wrapper ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree,
                             t8_eclass_t tree_class, t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme* scheme,
                             const int is_family, const int num_elements, t8_element_t* elements[])
{
  if (is_family && !scheme->elements_are_family (tree_class, elements)) {
    t8_global_errorf ("ERROR: The mesh handle does not support deleted elements.");
    return 0;  // No adaptation as default.
  }
  // Get static adapt context from the registry.
  // Via this, we can access the mesh handle and the user defined adapt callback that uses mesh handle functionality.
  auto* context = adapt_registry::get (forest_from);
  if (!context) {
    t8_global_errorf (
      "ERROR: Something went wrong while registering the adaptation callbacks. Please check your implementation.");
    return 0;  // No adaptation as default.
  }
  // Convert to index used in the mesh handle.
  const t8_locidx_t mesh_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  // Call the actual adapt callback stored in the context.
  return context->adapt_mesh (mesh_index, num_elements);
}

}  // namespace detail
}  // namespace t8_mesh_handle

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

/** \file interpolate.hxx
 * This file provides helper functionality to interpolate the element data of a \ref t8_mesh_handle::mesh
 * onto a newly adapted mesh according to a user defined callback.
 * During adaptation the element data attached to the old elements has to be transferred to the new (refined,
 * coarsened or unchanged) elements. t8code reports this old-to-new correspondence through
 * \ref t8_forest_iterate_replace, which expects a plain C callback of type \ref t8_forest_replace_t.
 * The helpers in this file bridge the gap between the mesh handle patterns and these functions,
 * following the same registry pattern as \ref adapt.hxx.
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <mesh_handle/mesh.hxx>
#include <memory>
#include <unordered_map>

namespace t8_mesh_handle
{

/** Namespace detail to hide implementation details from the user. */
namespace detail
{

/** Virtual base class for mesh interpolation contexts.
 * We need this base class and not only \ref mesh_interpolate_context for the \ref interpolate_registry.
 * interpolate_registry should not be templated because we need to access registered contexts in
 * \ref mesh_replace_callback_wrapper, where we do not know the type of the mesh. Therefore, we work with a map of
 * forests to instances of this (type erased) base class to remain template free.
 */
struct mesh_interpolate_context_base
{
  /** Virtual destructor for safe polymorphic deletion.
   */
  virtual ~mesh_interpolate_context_base () = default;

  /** Pure virtual callback to interpolate the element data of one set of old elements onto a set of new elements.
   * The indices refer to the flat, process local element numbering of the mesh handle (mesh handle id),
   * not to the tree local numbering used by the forest. The conversion is done in \ref mesh_replace_callback_wrapper.
   * \param [in] refine    -1 if a family in the old mesh got coarsened, 0 if the element was not touched,
   *                        1 if the element got refined.
   * \param [in] num_old    The number of outgoing (old) elements.
   * \param [in] first_old  The local mesh handle id of the first outgoing element in the old mesh.
   * \param [in] num_new    The number of incoming (new) elements.
   * \param [in] first_new  The local mesh handle id of the first incoming element in the new mesh.
   */
  virtual void
  interpolate (const int refine, const int num_old, const t8_locidx_t first_old, const int num_new,
               const t8_locidx_t first_new)
    = 0;
};

/** Templated mesh interpolation context holding the old and new mesh handle and the user defined callback.
 * Struct inherits from \ref mesh_interpolate_context_base and implements the virtual interpolate callback using the
 * two mesh handles and the callback.
 * Type erasure via the base class lets \ref interpolate_registry store this context without being templated on the
 * mesh type.
 * \tparam TMesh The mesh handle class.
 */
template <typename TMesh>
struct mesh_interpolate_context final: mesh_interpolate_context_base
{
  using callback_type =
    typename TMesh::internal_interpolate_callback_type; /**< The user defined interpolate callback type. */

  /** Constructor of the context with the old and new mesh handle and the user defined callback.
   * \param [in] mesh_old             The old mesh that is being adapted. Only read from during interpolation.
   * \param [in, out] mesh_new             The new mesh constructed from \a mesh_old. Written to during interpolation.
   * \param [in] interpolate_callback The interpolate callback. Moved into the context.
   */
  mesh_interpolate_context (const TMesh& mesh_old, TMesh& mesh_new, callback_type&& interpolate_callback)
    : m_mesh_old (mesh_old), m_mesh_new (mesh_new), m_callback (std::forward<callback_type> (interpolate_callback))
  {
  }

  /** Interpolation of one group of elements using the old and the new mesh and the user defined callback.
   * This function is called by \ref mesh_replace_callback_wrapper for each group.
   * \param [in] refine    -1 if a family got coarsened, 0 if the element was not touched, 1 if it got refined.
   * \param [in] num_old    The number of outgoing (old) elements.
   * \param [in] first_old  The local mesh handle id of the first outgoing element in the old mesh.
   * \param [in] num_new    The number of incoming (new) elements.
   * \param [in] first_new  The local mesh handle id of the first incoming element in the new mesh.
   */
  void
  interpolate (const int refine, const int num_old, const t8_locidx_t first_old, const int num_new,
               const t8_locidx_t first_new) override
  {
    // Check if the interpolate callback is set and call it using the old and new mesh handle.
    T8_ASSERTF (m_callback, "No interpolate callback set.");
    m_callback (m_mesh_old, m_mesh_new, refine, num_old, first_old, num_new, first_new);
  }

 private:
  const TMesh& m_mesh_old;        /**< The old mesh to read the element data from. */
  TMesh& m_mesh_new;              /**< The new mesh to write the interpolated element data to. */
  const callback_type m_callback; /**< The user defined interpolate callback. */
};

/** Registry pattern is used to register contexts, which provide access to the interpolate callback and the mesh
 * handles. This globally accessible static class is required to get the meshes and the callback in the forest
 * replace callback \ref mesh_replace_callback_wrapper, as the predefined \ref t8_forest_replace_t header does not
 * permit to pass these as function arguments. It uses the same idea as \ref adapt_registry.
 */
class interpolate_registry {
 public:
  /** Static function to register \a context using \a forest as identifier.
   * This makes the context publicly available through the registry.
   * \param [in] forest  The forest identifier. In our case, this is the old forest we interpolate data from.
   *                     ( Consistent to the adapt_registry.)
   * \param [in] context The context to register. Use unique pointer to ensure proper memory management and ownership.
   * \note We need the forest and not the mesh as key because this class must not be templated on the mesh type.
   */
  static void
  register_context (t8_forest_t forest, std::unique_ptr<mesh_interpolate_context_base> context)
  {
    auto& map = get_map ();
    auto [it, inserted] = map.emplace (forest, std::move (context));
    if (!inserted) {
      t8_global_errorf ("ERROR: Context already registered!");
    }
  }

  /** Static function to unregister a context using \a forest as identifier.
   * \param [in] forest The forest identifier. In our case, this is the old forest we interpolate from.
   */
  static void
  unregister_context (t8_forest_t forest)
  {
    auto& map = get_map ();
    [[maybe_unused]] const auto erased = map.erase (forest);
    T8_ASSERT (erased == 1);
  }

  /** Getter for a context using \a forest as identifier.
   * \param [in] forest The forest identifier. In our case, this is the old forest we interpolate from.
   * \return Pointer to the context registered with the id \a forest if found, nullptr otherwise.
   */
  static mesh_interpolate_context_base*
  get (t8_forest_t forest)
  {
    auto& map = get_map ();
    auto it = map.find (forest);
    return it != map.end () ? it->second.get () : nullptr;
  }

 private:
  /** Get the static map associating t8_forest_t with mesh_interpolate_context_base references.
   * We use a getter instead of private member variable to ensure single initialization.
   * \return Reference to the static unordered map of t8_forest_t to mesh_interpolate_context_base references.
   */
  static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_interpolate_context_base>>&
  get_map ()
  {
    static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_interpolate_context_base>> map;
    return map;
  }
};

/** Wrapper around the mesh handle interpolation functionality to be able to pass the callback to the classic replace
 * routine \ref t8_forest_iterate_replace of a forest. The function header fits the definition of
 * \ref t8_forest_replace_t. 
 * \param [in] forest_old      The forest that is adapted fitting the old element data.
 * \param [in, out] forest_new The forest that is newly constructed from \a forest_old. 
 *                             Data will be interpolated to match this forest.
 * \param [in] which_tree      The local tree containing \a first_outgoing and \a first_incoming.
 * \param [in] tree_class      Unused; The eclass of the local tree containing \a first_outgoing and \a first_incoming.
 * \param [in] scheme          Unused; The scheme of the forest.
 * \param [in] refine          -1 if a family got coarsened, 0 if the element was not touched, 1 if it got refined.
 * \param [in] num_outgoing    The number of outgoing elements.
 * \param [in] first_outgoing  The tree local index of the first outgoing element.
 *                             0 <= first_outgoing < which_tree->num_elements
 * \param [in] num_incoming    The number of incoming elements.
 * \param [in] first_incoming  The tree local index of the first incoming element.
 *                             0 <= first_incom < new_which_tree->num_elements
 */
inline void
mesh_replace_callback_wrapper (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                               [[maybe_unused]] const t8_eclass_t tree_class,
                               [[maybe_unused]] const t8_scheme_c* scheme, const int refine, const int num_outgoing,
                               const t8_locidx_t first_outgoing, const int num_incoming,
                               const t8_locidx_t first_incoming)
{
  // Get the static interpolate context from the registry.
  // Via this, we can access the old and new mesh handle and the user defined interpolate callback.
  auto* context = interpolate_registry::get (forest_old);
  if (!context) {
    t8_global_productionf ("Interpolate context not found. Did you forget to register it?");
    return;
  }

  // Convert the tree local indices reported by the forest to the flat, process local indices used in the mesh handle
  // (the mesh handle id).
  const t8_locidx_t first_old_global = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_outgoing;
  const t8_locidx_t first_new_global = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;
  // Call the actual interpolate callback stored in the context.
  context->interpolate (refine, num_outgoing, first_old_global, num_incoming, first_new_global);
}

}  // namespace detail
}  // namespace t8_mesh_handle

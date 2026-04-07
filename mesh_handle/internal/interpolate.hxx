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
 * This file provides helper functionality to interpolate element data for a \ref t8_mesh_handle::mesh
 * according to a user defined callback.
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <mesh_handle/mesh.hxx>
#include <memory>
#include <unordered_map>

namespace t8_mesh_handle
{
namespace detail
{

/**
 * Base class for interpolation contexts (type-erased).
 */
struct mesh_interpolate_context_base
{
  virtual ~mesh_interpolate_context_base () = default;

  /**
   * Virtual interpolation callback.
   */
  virtual void
  interpolate (const int refine, const int num_old, const t8_locidx_t first_old, const int num_new,
               const t8_locidx_t first_new)
    = 0;
};

/**
 * Concrete interpolation context storing mesh handles + user callback.
 */
template <typename TMesh>
struct mesh_interpolate_context final: mesh_interpolate_context_base
{
  using callback_type = typename TMesh::interpolate_callback_type;

  mesh_interpolate_context (const TMesh& mesh_old, TMesh& mesh_new, callback_type interpolate_callback)
    : m_mesh_old (mesh_old), m_mesh_new (mesh_new), m_callback (std::move (interpolate_callback))
  {
  }

  void
  interpolate (const int refine, const int num_old, const t8_locidx_t first_old, const int num_new,
               const t8_locidx_t first_new) override
  {
    T8_ASSERTF (m_callback, "No interpolate callback set.");
    m_callback (m_mesh_old, m_mesh_new, refine, num_old, first_old, num_new, first_new);
  }

 private:
  const TMesh& m_mesh_old;
  TMesh& m_mesh_new;
  const callback_type m_callback;
};

/**
 * Registry for interpolate contexts (same idea as adapt_registry).
 */
class interpolate_registry {
 public:
  /** Need forest and not mesh for key because this cannot have a template parameter. */
  static void
  register_context (t8_forest_t forest, std::unique_ptr<mesh_interpolate_context_base> context)
  {
    auto& map = get_map ();
    auto [it, inserted] = map.emplace (forest, std::move (context));
    if (!inserted) {
      t8_global_productionf ("Interpolate context already registered!");
    }
  }

  static void
  unregister_context (t8_forest_t forest)
  {
    auto& map = get_map ();
    [[maybe_unused]] const auto erased = map.erase (forest);
    T8_ASSERT (erased == 1);
  }

  static mesh_interpolate_context_base*
  get (t8_forest_t forest)
  {
    auto& map = get_map ();
    auto it = map.find (forest);
    return it != map.end () ? it->second.get () : nullptr;
  }

 private:
  static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_interpolate_context_base>>&
  get_map ()
  {
    static std::unordered_map<t8_forest_t, std::unique_ptr<mesh_interpolate_context_base>> map;
    return map;
  }
};

/**
 * Wrapper matching t8_forest_replace_t
 */
inline void
mesh_replace_callback_wrapper (t8_forest_t forest_old, [[maybe_unused]] t8_forest_t forest_new, t8_locidx_t which_tree,
                               const t8_eclass_t tree_class, const t8_scheme_c* scheme, const int refine,
                               const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
                               const t8_locidx_t first_incoming)
{
  auto* context = interpolate_registry::get (forest_old);
  if (!context) {
    t8_global_productionf ("Interpolate context not found. Did you forget to register it?");
    return;
  }

  // Convert tree-local indices to mesh-global indices
  const t8_locidx_t first_old_global = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_outgoing;
  const t8_locidx_t first_new_global = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;
  context->interpolate (refine, num_outgoing, first_old_global, num_incoming, first_new_global);
}

}  // namespace detail
}  // namespace t8_mesh_handle

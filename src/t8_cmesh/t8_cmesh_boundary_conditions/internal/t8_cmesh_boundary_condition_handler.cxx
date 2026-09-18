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

/**
 * \file t8_cmesh_boundary_condition_handler.cxx
 * Implementation context of \ref t8_cmesh_boundary_condition_handler.hxx
 */

#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler.hxx>
#include <t8_cmesh/t8_cmesh_boundary_conditions/internal/t8_cmesh_boundary_condition_handler_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_types.h>
#include <t8_cmesh/t8_cmesh_internal/t8_cmesh_stash.h>

using namespace detail;

#if T8_ENABLE_DEBUG
int
t8_cmesh_boundary_condition_handler::verify () const
{
  T8_ASSERT (!t8_cmesh_is_committed (m_cmesh, 0));

  /* Retrieve eclasses and boundary conditions from stash */
  std::vector<std::pair<t8_gloidx_t, std::span<boundary_condition_hash>>> boundary_conditions
    = t8_stash_extract_attribute_list<boundary_condition_hash> (m_cmesh->stash, t8_get_package_id (),
                                                                T8_CMESH_BOUNDARY_CONDITION_ATTRIBUTE_KEY);
  std::vector<std::pair<t8_gloidx_t, t8_eclass>> eclasses = t8_stash_extract_eclasses (m_cmesh->stash);

  /* Check that every tree has a bc */
  if (eclasses.size () != boundary_conditions.size ()) {
    t8_errorf ("ERROR: Number of cmesh eclasses does not match the number of cmesh boundary conditions.\n"
               "If boundary conditions are applied, all trees have to get them assigned. They can be left empty.\n");
    return false;
  }

  /* Sort the bcs and eclasses by tree id. */
  std::sort (boundary_conditions.begin (), boundary_conditions.end (),
             [] (const auto &a, const auto &b) { return a.first < b.first; });
  std::sort (eclasses.begin (), eclasses.end (), [] (const auto &a, const auto &b) { return a.first < b.first; });

  auto eclasses_it = eclasses.cbegin ();
  auto bc_it = boundary_conditions.cbegin ();
  for (; eclasses_it != eclasses.end (); ++eclasses_it, ++bc_it) {
    /* Make sure, that the eclasses and bcs tree id are the same.  */
    if (eclasses_it->first != bc_it->first) {
      t8_errorf ("ERROR: eclass tree id and boundary condition tree id do not match.\nProbably the boundary condition "
                 "got a wrong tree id. Boundary condition tree id: %li \n",
                 bc_it->first);
      return false;
    }
    /* Make sure that the amount of registered bcs matches th number of tree faces. */
    if (static_cast<size_t> (t8_eclass_num_faces[eclasses_it->second]) != bc_it->second.size ()) {
      t8_errorf ("ERROR: Tree %li has a mismatch in its face count and boundary condition count.\n"
                 "Face count: %i, boundary condition count: %li\n",
                 eclasses_it->first, t8_eclass_num_faces[eclasses_it->second], bc_it->second.size ());
      return false;
    }
    /* Check, that the hash of the bc attribute is present in the bc map. This failure cannot be caused by the user. */
    for (const auto bc : bc_it->second) {
      if (!get_boundary_condition_name_safe (bc).has_value ()) {
        t8_errorf ("ERROR: Could not find name of applied boundary condition. You probably found a bug.\n");
        return false;
      }
    }
  }
  return true;
}

#endif /* T8_ENABLE_DEBUG */

int
t8_cmesh_boundary_condition_handler::get_boundary_condition_attribute_key () const
{
  return T8_CMESH_BOUNDARY_CONDITION_ATTRIBUTE_KEY;
}

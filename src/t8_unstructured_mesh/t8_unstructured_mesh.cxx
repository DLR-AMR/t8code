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

/** \file t8_unstructured_mesh.cxx
 *  This file implements functions defined in the unstructured mesh element class.
 */

#include <t8_unstructured_mesh/t8_unstructured_mesh.hxx>

int
t8_unstructured_mesh_element::get_level ()
{
  if (m_unstructured_mesh->m_level_cache.empty ()) {
    const t8_scheme* scheme = t8_forest_get_scheme (m_unstructured_mesh->m_forest);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (m_unstructured_mesh->m_forest, m_tree_id);
    const t8_element_t* element
      = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
    return scheme->element_get_level (tree_class, element);
  }
  else {
    return m_unstructured_mesh->m_level_cache[m_tree_id][m_element_id];
  }
}

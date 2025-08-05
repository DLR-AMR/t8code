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

#include <t8_unstructured_mesh/t8_unstructured_mesh.hxx>

t8_unstructured_mesh::Element_Iterator::Element_Iterator (t8_unstructured_mesh* unstructured_mesh,
                                                          t8_locidx_t current_tree_id, t8_locidx_t current_element_id)
  : m_current_tree_id (current_tree_id), m_current_element_id (current_element_id),
    m_unstructured_mesh (unstructured_mesh)
{
  m_num_local_trees = t8_forest_get_num_local_trees (m_unstructured_mesh->m_forest);
  if (m_num_local_trees > m_current_tree_id) {
    m_num_elements_current_tree
      = t8_forest_get_tree_num_leaf_elements (m_unstructured_mesh->m_forest, m_current_tree_id);
  }
  else {
    m_num_elements_current_tree = 0;
  }
}

t8_unstructured_mesh::Element_Iterator::reference
t8_unstructured_mesh::Element_Iterator::operator* () const
{
  return *new t8_unstructured_mesh_element (m_unstructured_mesh, m_current_tree_id, m_current_element_id);
}

t8_unstructured_mesh::Element_Iterator::pointer
t8_unstructured_mesh::Element_Iterator::operator->() const
{
  return new t8_unstructured_mesh_element (m_unstructured_mesh, m_current_tree_id, m_current_element_id);
}

t8_unstructured_mesh::Element_Iterator&
t8_unstructured_mesh::Element_Iterator::operator++ ()
{
  if (m_current_element_id < m_num_elements_current_tree - 1) {
    m_current_element_id++;
  }
  else {
    m_current_element_id = 0;
    m_current_tree_id++;
    if (m_num_local_trees > m_current_tree_id) {
      m_num_elements_current_tree
        = t8_forest_get_tree_num_leaf_elements (m_unstructured_mesh->m_forest, m_current_tree_id);
    }
    else {
      m_num_elements_current_tree = 0;
    }
  }
  return *this;
}

t8_unstructured_mesh::Element_Iterator
t8_unstructured_mesh::Element_Iterator::operator++ (int)
{
  t8_unstructured_mesh::Element_Iterator tmp = *this;
  ++(*this);
  return tmp;
}

bool
t8_unstructured_mesh::Element_Iterator::operator== (const t8_unstructured_mesh::Element_Iterator& other_iterator) const
{
  return m_unstructured_mesh->m_forest == other_iterator.m_unstructured_mesh->m_forest
         && m_current_tree_id == other_iterator.m_current_tree_id
         && m_current_element_id == other_iterator.m_current_element_id;
}
// Not needed in C++20 but for completion.
bool
t8_unstructured_mesh::Element_Iterator::operator!= (const t8_unstructured_mesh::Element_Iterator& other_iterator) const
{
  return !(*this == other_iterator);
}

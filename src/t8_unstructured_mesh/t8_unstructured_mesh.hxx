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
#ifndef T8_UNSTRUCTURED_MESH_HXX
#define T8_UNSTRUCTURED_MESH_HXX
#include <t8.h>
#include <t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <iterator>
#include <cstddef>

/** TODO: Maybe templated class with template unstructured mesh element? Then we could also use a child of 
unstructured_elements here and construct their type. We want to derive from unstructured element if we want to add 
a further property of the element that can be accessed. */
template <class unstructured_mesh_element>
class t8_unstructured_mesh {

  friend unstructured_mesh_element;

 public:
  /** TODO*/
  t8_unstructured_mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
  }

  /** \brief This iterator should iterate over all (local) elements.
  */
  struct Element_Iterator
  {
    /* Design choice: This Iterator is part of the unstructured mesh as it is */
    using iterator_category = std::forward_iterator_tag;  //TODO: do we maybe need a bi-directional Iterator?
    using difference_type = std::ptrdiff_t;
    using value_type = unstructured_mesh_element;
    using pointer = value_type*;
    using reference = value_type&;
    // Constructor.
    Element_Iterator (t8_unstructured_mesh* unstructured_mesh, t8_locidx_t current_tree_id,
                      t8_locidx_t current_element_id)
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

    reference
    operator* () const
    {
      return *new unstructured_mesh_element (m_unstructured_mesh, m_current_tree_id, m_current_element_id);
    }
    pointer
    operator->() const
    {
      return new unstructured_mesh_element (m_unstructured_mesh, m_current_tree_id, m_current_element_id);
    }

    Element_Iterator&
    operator++ ()
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

    Element_Iterator
    operator++ (int)
    {
      t8_unstructured_mesh::Element_Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool
    operator== (const Element_Iterator& other_iterator) const
    {
      return m_unstructured_mesh->m_forest == other_iterator.m_unstructured_mesh->m_forest
             && m_current_tree_id == other_iterator.m_current_tree_id
             && m_current_element_id == other_iterator.m_current_element_id;
    }

    // Not needed in C++20 but for completion.
    bool
    operator!= (const Element_Iterator& other_iterator) const
    {
      return !(*this == other_iterator);
    }

   private:
    t8_locidx_t m_current_tree_id, m_current_element_id;
    t8_unstructured_mesh* m_unstructured_mesh;
    t8_locidx_t m_num_local_trees, m_num_elements_current_tree;
    unstructured_mesh_element* current_element;
  };

  /**TODO*/
  inline Element_Iterator
  begin ()
  {
    return Element_Iterator (this, 0, 0);
  }

  /**
 * TODO
 */
  inline Element_Iterator
  end ()
  {
    return Element_Iterator (this, t8_forest_get_num_local_trees (m_forest), 0);
  }

 private:
  t8_forest_t m_forest;
};

class t8_unstructured_mesh_element {
 public:
  t8_unstructured_mesh_element (t8_unstructured_mesh<t8_unstructured_mesh_element>* unstruct, t8_locidx_t tree_id,
                                t8_locidx_t element_id)
    : m_tree_id (tree_id), m_element_id (element_id), m_unstructured_mesh (unstruct)
  {
  }

  int
  get_level ();

 private:
  t8_locidx_t m_tree_id, m_element_id;
  t8_unstructured_mesh<t8_unstructured_mesh_element>* m_unstructured_mesh;
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */

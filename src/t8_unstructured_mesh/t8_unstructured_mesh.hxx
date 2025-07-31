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

//TODO: Inspiration by t8_element_array_iterator
class t8_unstructured_mesh {
 public:
  t8_unstructured_mesh (t8_forest_t input_forest): forest (input_forest)
  {
  }

  /** \brief This iterator should iterate over all (local) elements.
  */
  struct Element_Iterator
  {
    using iterator_category = std::forward_iterator_tag;  //TODO: do we maybe need a bidrirecIterator?
    using difference_type = std::ptrdiff_t;
    using value_type = t8_element_t*;
    using pointer = value_type*;
    using reference = value_type&;
    // Constructor.
    Element_Iterator (t8_forest_t forest, t8_locidx_t current_tree_id, t8_locidx_t current_element_id)
      : m_current_tree_id (current_tree_id), m_current_element_id (current_element_id), m_forest (forest)
    {
      m_num_local_trees = t8_forest_get_num_local_trees (m_forest);
      m_num_elements_current_tree = t8_forest_get_tree_num_leaf_elements (m_forest, m_current_tree_id);
      m_scheme = t8_forest_get_scheme (forest);
    }

    const t8_element_t*
    operator* () const
    {
      const t8_element_t* elem = t8_forest_get_leaf_element_in_tree (m_forest, m_current_tree_id, m_current_element_id);
      if (elem == nullptr) {
        SC_ABORT ("not implemented yet");
      }
      return elem;
    }

    // Prefix version of ++.
    Element_Iterator&
    operator++ ()
    {
      if (m_current_element_id < m_num_elements_current_tree - 1) {
        m_current_element_id++;
      }
      else {
        m_current_element_id = 0;
        m_current_tree_id++;
        m_num_elements_current_tree
          = t8_forest_get_tree_num_leaf_elements (m_forest, m_current_tree_id);  //Problem for last element, right?
      }
      return *this;
    }

    // Postfix version of ++.
    Element_Iterator
    operator++ (int)
    {
      Element_Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool
    operator== (const Element_Iterator& other_iterator) const
    {
      return m_forest == other_iterator.m_forest && m_current_tree_id == other_iterator.m_current_tree_id
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
    t8_forest_t m_forest;
    t8_locidx_t m_num_local_trees, m_num_elements_current_tree;
    const t8_scheme* m_scheme;
  };

  /**TODO*/
  inline Element_Iterator
  element_begin ()
  {
    return Element_Iterator (forest, 0, 0);
  }

  /**
 * TODO
 */
  inline Element_Iterator
  element_end ()
  {
    return Element_Iterator (forest, t8_forest_get_num_local_trees (forest) - 1,
                             t8_forest_get_tree_num_leaf_elements (forest, t8_forest_get_num_local_trees (forest) - 1));
  }

 private:
  t8_forest_t forest;
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */

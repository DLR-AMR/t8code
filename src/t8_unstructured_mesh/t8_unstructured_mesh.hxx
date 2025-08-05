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

class t8_unstructured_mesh_element;
/** TODO*/
class t8_unstructured_mesh {

  friend class t8_unstructured_mesh_element;

 public:
  /** TODO*/
  t8_unstructured_mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
  }

  /** \brief This iterator should iterate over all (local) elements.
  */
  struct Element_Iterator
  {
    using iterator_category = std::forward_iterator_tag;  //TODO: do we maybe need a bidrirecIterator?
    using difference_type = std::ptrdiff_t;
    using value_type = t8_unstructured_mesh_element;
    using pointer = value_type*;
    using reference = value_type&;
    // Constructor.
    Element_Iterator (t8_unstructured_mesh* unstructured_mesh, t8_locidx_t current_tree_id,
                      t8_locidx_t current_element_id);

    reference
    operator* () const;
    pointer
    operator->() const;

    Element_Iterator&
    operator++ ();
    Element_Iterator
    operator++ (int);

    bool
    operator== (const Element_Iterator& other) const;
    bool
    operator!= (const Element_Iterator& other) const;

   private:
    t8_locidx_t m_current_tree_id, m_current_element_id;
    t8_unstructured_mesh* m_unstructured_mesh;
    t8_locidx_t m_num_local_trees, m_num_elements_current_tree;
    t8_unstructured_mesh_element* current_element;
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
  t8_unstructured_mesh_element (t8_unstructured_mesh* unstruct, t8_locidx_t tree_id, t8_locidx_t element_id)
    : m_tree_id (tree_id), m_element_id (element_id), m_unstructured_mesh (unstruct)
  {
  }

  int
  get_level ()
  {
    const t8_scheme* scheme = t8_forest_get_scheme (m_unstructured_mesh->m_forest);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (m_unstructured_mesh->m_forest, m_tree_id);
    const t8_element_t* element
      = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
    return scheme->element_get_level (tree_class, element);
  }

 private:
  t8_locidx_t m_tree_id, m_element_id;
  t8_unstructured_mesh* m_unstructured_mesh;
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */

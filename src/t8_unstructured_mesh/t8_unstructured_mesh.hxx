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

/** \file t8_unstructured_mesh.hxx
 * Definition of the unstructured mesh class and related functionality.
 */

#ifndef T8_UNSTRUCTURED_MESH_HXX
#define T8_UNSTRUCTURED_MESH_HXX
#include <t8.h>

#include <t8_element.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <iterator>
#include <memory>
#include <vector>

/* Forward declaration of the default unstructured element class because used as default template parameter in
 * t8_unstructured_mesh but t8_unstructured_mesh is also needed to define the element class.
 */
class t8_unstructured_mesh_element;

/**
 * Wrapper for a forest to be handled like an unstructured mesh object. 
 * The default class t8_unstructured_mesh_element provides access to the default functionality needed.
 * In the unstructured mesh class, you can decide which parameters should be cached and which should be calculated on the fly. 
 * Per default, the parameters will be calculated and not cached, please call the related functions to cache variables.
 * If you want to access more element parameters than the default ones, that should not be cached, you can write a derived class of 
 * t8_unstructured_mesh_element and use the derived class as template parameter. If the additional variable(s) should be cached,
 *  you may also write a derived class of t8_unstructured_mesh.
 * \tparam unstructured_mesh_element The element class that should be used for the unstructured mesh elements. 
 */
template <class TUnstructuredMeshElement = t8_unstructured_mesh_element>
class t8_unstructured_mesh {
  using element_vector = std::vector<TUnstructuredMeshElement>;

 public:
  // Declare unstructured mesh element as friend such that the forest and cached variables can be accessed.
  friend TUnstructuredMeshElement;

  /** 
   * Constructor for an unstructured mesh. 
   * \param [in] input_forest The forest from which the unstructured mesh should be created. 
   */
  t8_unstructured_mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
    m_scheme = t8_forest_get_scheme (m_forest);
    update_elements ();
  }

  void
  update_elements ()
  {
    if (!m_elements.empty ()) {
      m_elements.clear ();
    }
    m_num_local_trees = t8_forest_get_num_local_trees (m_forest);
    for (t8_locidx_t itree = 0; itree < m_num_local_trees; ++itree) {
      const t8_locidx_t num_elems = t8_forest_get_tree_num_leaf_elements (m_forest, itree);
      element_vector temp;
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        temp.push_back (TUnstructuredMeshElement (this, itree, ielem));
      }
      m_elements.push_back (temp);
    }
  }

  /** This forward iterator iterates over all (local) elements of the unstructured mesh.
   */
  struct t8_unstructured_iterator
  {
    /* Design choice: This Iterator is part of the unstructured mesh as it is strongly connected to the unstructured mesh 
     * and we should not need derived iterators. 
     */
    using iterator_category = std::forward_iterator_tag;  //TODO: do we maybe need a bi-directional Iterator?
    using difference_type = std::ptrdiff_t;
    using value_type = TUnstructuredMeshElement;
    using pointer = value_type*;
    using reference = value_type&;

    /** 
     * Constructor for the element iterator. 
     * \param [in] unstructured_mesh Pointer to the unstructured mesh the iterator should be created for. 
     * \param [in] current_tree_id Initial tree id of the iterator. 
     * \param [in] current_element_id Initial element id in the tree of the iterator. 
     */
    t8_unstructured_iterator (t8_unstructured_mesh* unstructured_mesh, t8_locidx_t current_tree_id,
                              t8_locidx_t current_element_id)
      : m_unstructured_mesh (unstructured_mesh)
    {
      m_outer_iterator = m_unstructured_mesh->m_elements.begin () + current_tree_id;
      m_inner_iterator = m_outer_iterator->begin () + current_element_id;
    }

    /**
     * Dereference the iterator to access the unstructured mesh element.
     * \return Reference to the current unstructured mesh element.
     */
    reference
    operator* () const
    { /* Define new unstructured mesh element instead of caching the current one because it is possible that the iterator 
       * should point to an element following the last (local) element, e.g., end().
       */
      return (*m_inner_iterator);
    }

    /**
     * Access member of the current unstructured mesh element.
     * \return Pointer to the current unstructured mesh element.
     */
    pointer
    operator->() const
    {
      return &(*m_inner_iterator);
    }

    /**
     * Prefix-increment the iterator to point to the next element.
     * \return Reference to the incremented iterator.
     */
    t8_unstructured_iterator&
    operator++ ()
    {
      if (m_inner_iterator == (*m_outer_iterator).end () - 1) {
        m_outer_iterator++;
        m_inner_iterator = m_outer_iterator->begin ();
      }
      else {
        m_inner_iterator++;
      }
      return *this;
    }

    /**
     * Post-increment the iterator.
     * \return Iterator before increment.
     */
    t8_unstructured_iterator
    operator++ (int)
    {
      t8_unstructured_mesh::t8_unstructured_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    /**
     * Equality comparison.
     * 
     * \param [in] other_iterator Another iterator to compare.
     * \return True if both iterators point to the same element, false otherwise.
     */
    bool
    operator== (const t8_unstructured_iterator& other_iterator) const
    {
      return m_unstructured_mesh->m_forest == other_iterator.m_unstructured_mesh->m_forest
             && m_outer_iterator == other_iterator.m_outer_iterator
             && m_inner_iterator == other_iterator.m_inner_iterator;
    }

    /**
     * Inequality comparison operator.
     * This operator is not needed in C++20 but for completion.
     * 
     * \param [in] other_iterator Another iterator to compare.
     * \return True if both iterators point to different elements, false otherwise.
     */
    bool
    operator!= (const t8_unstructured_iterator& other_iterator) const
    {
      return !(*this == other_iterator);
    }

   private:
    t8_unstructured_mesh* m_unstructured_mesh; /**< The unstructured mesh the iterator is defined for. */
    std::vector<element_vector>::iterator m_outer_iterator;
    element_vector::iterator m_inner_iterator;
  };

  /**
   * Returns an iterator to the first (local) unstructured mesh element.
   */
  inline t8_unstructured_iterator
  begin ()
  {
    return t8_unstructured_iterator (this, 0, 0);
  }

  /**
   * Returns an iterator to an unstructured mesh element following the last (local) element of the unstructured mesh.
   */
  inline t8_unstructured_iterator
  end ()
  {
    return t8_unstructured_iterator (this, t8_forest_get_num_local_trees (m_forest), 0);
  }

 private:
  t8_forest_t m_forest; /**< The forest the unstructured mesh should be defined for. */
  const t8_scheme* m_scheme;
  t8_locidx_t m_num_local_trees;
  std::vector<element_vector> m_elements;
};

/** 
 * Default element of an unstructured mesh. 
 * For this element, the following properties can be accessed: Level, TODO.
 */
class t8_unstructured_mesh_element {
  /* Design choice: Decided to not define the class inside of \ref t8_unstructured_mesh although the classes are strongly connected,
* because the class also would not have access to private members and inheritance of the element class would be complicated.
*/
 public:
  t8_unstructured_mesh_element (t8_unstructured_mesh<t8_unstructured_mesh_element>* unstructured_mesh,
                                t8_locidx_t tree_id, t8_locidx_t element_id)
    : m_tree_id (tree_id), m_element_id (element_id), m_unstructured_mesh (unstructured_mesh)
  {
  }

  /**
   * Getter for the refinement level of the unstructured mesh element.
   * \return Refinement level of the unstructured mesh element.
   */
  t8_element_level
  get_level ();

 private:
  t8_locidx_t m_tree_id,
    m_element_id; /**< The tree id and the element id of the element in the forest defined in the unstructured mesh. */
  t8_unstructured_mesh<t8_unstructured_mesh_element>*
    m_unstructured_mesh; /**< Pointer to the unstructured mesh the element is defined for. */
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */

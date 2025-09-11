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
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <iterator>
#include <memory>
#include <vector>
#include <tuple>
#include <t8_unstructured_mesh/t8_element_competences.hxx>

/* Forward declaration of the default unstructured element class used as default template parameter in
 * t8_unstructured_mesh.
 */
template <template <typename> class... TCompetence>
class t8_unstructured_mesh_element;

/**
 * Wrapper for a forest that enables it to be handled like an unstructured mesh object.
 * \tparam unstructured_mesh_element The element class that should be used for the unstructured mesh elements. 
 * This template parameter defines which element functionality is available and if it is cached or calculated.
 */
template <class TUnstructuredMeshElement = t8_unstructured_mesh_element<>>
class t8_unstructured_mesh {
  using element_vector = std::vector<TUnstructuredMeshElement>;

 public:
  // Declare unstructured mesh element as friend such that the forest can be accessed.
  friend TUnstructuredMeshElement;

  /** 
   * Constructor for an unstructured mesh. 
   * \param [in] input_forest The forest from which the unstructured mesh should be created. 
   */
  t8_unstructured_mesh (t8_forest_t input_forest): m_forest (input_forest)
  {
    update_elements ();
  }

  /** 
   * Update the storage of the unstructured mesh elements according to the current forest. 
   * Can be used for example after the forest is adapted.  
   */
  void
  update_elements ()
  {
    // Clear the element vector if already created.
    if (!m_elements.empty ()) {
      m_elements.clear ();
    }
    // Iterate through forest elements and fill the element vector with newly created unstructured mesh elements.
    for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (m_forest); ++itree) {
      const t8_locidx_t num_elems = t8_forest_get_tree_num_leaf_elements (m_forest, itree);
      element_vector temp;
      temp.reserve (num_elems);
      for (t8_locidx_t ielem = 0; ielem < num_elems; ++ielem) {
        temp.emplace_back (this, itree, ielem);
      }
      m_elements.push_back (std::move (temp));
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
      // Check if the outer iterator is pointing to an valid vector.
      if (m_outer_iterator != m_unstructured_mesh->m_elements.end ()) {
        m_inner_iterator = m_outer_iterator->begin () + current_element_id;
      }
      else {
        // If the outer iterator points to the end of the vector, define the current position of the inner
        // iterator to end() of the last vector in the element vector. This is also the natural way for increment.
        m_inner_iterator = (m_unstructured_mesh->m_elements.end () - 1)->end ();
      }
    }

    /**
     * Dereference the iterator to access the unstructured mesh element.
     * \return Reference to the current unstructured mesh element.
     */
    reference
    operator* () const
    {
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
      // Check if the iterator is at the end of the current tree vector. Next vector if yes and inner iterator increment if not.
      if (m_inner_iterator == (*m_outer_iterator).end () - 1) {
        m_outer_iterator++;
        // If the outer iterator does not point to a valid vector anymore,
        // set the inner iterator also to end(), else, to begin() of the new tree vector.
        if (m_outer_iterator != m_unstructured_mesh->m_elements.end ()) {
          m_inner_iterator = m_outer_iterator->begin ();
        }
        else {
          m_inner_iterator++;
        }
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
    typename std::vector<element_vector>::iterator
      m_outer_iterator;                                 /**< The iterator for the outer vector of the element vector. */
    typename element_vector::iterator m_inner_iterator; /**< The iterator for the inner vector of the element vector. */
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

  /**
   * Getter for the forest the unstructured mesh is defined for.
   */
  t8_forest_t
  get_forest ()
  {
    return m_forest;
  }

 private:
  t8_forest_t m_forest; /**< The forest the unstructured mesh should be defined for. */
  std::vector<element_vector>
    m_elements; /**< Vector storing the unstructured mesh elements. One element vector per (local) tree. */
};

/** 
 * Unstructured mesh element class. 
 * The unstructured element without specified template parameters provides default implementations for basic functionality 
 * as accessing the refinement level or the centroid. With this implementation, the functionality is calculated each time
 * the function is called. 
 * Use the competences defined in t8_element_competences.hxx as template parameter to cache the functionality instead of 
 * calculating them each time. 
 * To add functionality to the element, you can simply write you own competence class and give it as a template parameter.
 * You can access the functions implemented in your competence via the element. 
 *
 * The inheritance pattern is inspired by the \ref T8Type class.
 * We decided to use this structure 1.) to be able to add new functionality easily and 
 *    2.) for the cached options to keep the number of class member variables of the default to a minimum to safe memory.
 * The choice between calculate and cache is a tradeoff between runtime and memory usage. 
 *
 * \tparam The competences you want to add to the default functionality of the element.
 */
template <template <typename> class... TCompetence>
class t8_unstructured_mesh_element: public TCompetence<t8_unstructured_mesh_element<TCompetence...>>... {
  using SelfType = t8_unstructured_mesh_element<TCompetence...>;

  // --- Variables to check which functionality is defined in TCompetence. ---
  // Checks if one of the competences (like CacheLevel) defines the function get_level_cached().
  static constexpr bool get_level_defined
    = (false || ... || requires (TCompetence<SelfType>& competence) { competence.get_level_cached (); });
  static constexpr bool get_centroid_defined
    = (false || ... || requires (TCompetence<SelfType>& competence) { competence.get_centroid_cached (); });

 public:
  /**
   * Constructor of the unstructured mesh element.
   * \param unstructured_mesh Reference to the unstructured mesh the element should belong to.
   * \param tree_id The tree id of the element in the forest defining the unstructured mesh.
   * \param element_id The element id of the element in the forest defining the unstructured mesh.
   */
  t8_unstructured_mesh_element (t8_unstructured_mesh<SelfType>* unstructured_mesh, t8_locidx_t tree_id,
                                t8_locidx_t element_id)
    : m_tree_id (tree_id), m_element_id (element_id), m_unstructured_mesh (unstructured_mesh)
  {
  }

  // --- Functionality of the element. In each function, it is checked if a cached version exists (and is used then). ---

  /**
   * Getter for the refinement level of the unstructured mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates the refinement level if not.
   * \return Refinement level of the unstructured mesh element.
   */
  t8_element_level
  get_level ()
  {
    if constexpr (get_level_defined) {
      return this->get_level_cached ();
    }
    else {
      const t8_eclass_t tree_class = t8_forest_get_tree_class (m_unstructured_mesh->m_forest, m_tree_id);
      const t8_element_t* element
        = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
      return t8_forest_get_scheme (m_unstructured_mesh->m_forest)->element_get_level (tree_class, element);
    }
  }

  /**
   * Getter for the center of mass of the unstructured mesh element.
   * This function uses the cached version defined in TCompetence if available and calculates if not.
   * \return Coordinates of the center.
   */
  double*
  get_centroid ()
  {
    if constexpr (get_centroid_defined) {
      return this->get_centroid_cached ();
    }
    else {
      double* coordinates = new double[t8_forest_get_dimension (m_unstructured_mesh->m_forest)];
      const t8_element_t* element
        = t8_forest_get_leaf_element_in_tree (m_unstructured_mesh->m_forest, m_tree_id, m_element_id);
      t8_forest_element_centroid (m_unstructured_mesh->m_forest, m_tree_id, element, coordinates);
      return coordinates;
    }
  }

  //--- Getter for the member variables. ---
  /**
   * Getter for the tree id of the unstructured mesh element.
   */
  t8_locidx_t
  get_tree_id ()
  {
    return m_tree_id;
  }

  /**
   * Getter for the element id of the unstructured mesh element.
   */
  t8_locidx_t
  get_element_id ()
  {
    return m_element_id;
  }

  /**
   * Getter for the unstructured mesh to which the unstructured mesh element is belonging.
   * \return Reference to the unstructured mesh.
   */
  t8_unstructured_mesh<SelfType>*
  get_unstructured_mesh ()
  {
    return m_unstructured_mesh;
  }

 private:
  t8_locidx_t m_tree_id,
    m_element_id; /**< The tree id and the element id of the element in the forest defined in the unstructured mesh. */
  t8_unstructured_mesh<SelfType>*
    m_unstructured_mesh; /**< Pointer to the unstructured mesh the element is defined for. */
};

#endif /* !T8_UNSTRUCTURED_MESH_HXX */

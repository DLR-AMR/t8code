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

/** \file t8_forest_adapt.hxx 
 * Definition of a C++ class for the adaptation routine to provide a flexible way of implementing
 * different kinds of adaptation strategies.
 */

#ifndef T8_FOREST_ADAPT_HXX
#define T8_FOREST_ADAPT_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_types.h>

#include <functional>

/**
 * Namespace for adaptation related classes and functions.
 */

 /* TODO rename to t8_forest_adapt as soon as it is not used as function name anymore.  */
namespace t8_forest_adapt_namespace {
  /** The action to be taken on an element during adaptation. 
   * COARSEN: The element should be coarsened.
   * KEEP: The element should remain as is.
   * REFINE: The element should be refined.
  */
  enum adapt_action {
    COARSEN = -1,
    KEEP = 0,
    REFINE = 1
  };

  /**
   * Callback function type for element adaptation.
   */
  using element_callback = std::function<adapt_action(const t8_forest_t forest, const t8_locidx_t ltreeid,
                                                     const t8_element_t *element, const t8_scheme *scheme, const t8_eclass_t tree_class)>;

  /**  * Function to manipulate elements based on the specified adaptation action.
   * \tparam action The adaptation action to be performed.
   * \param [in, out] elements         The element array to be modified.
   * \param [in]     elements_from    The source element array.
   * \param [in]     scheme           The element scheme.
   * \param [in]     tree_class       The eclass of the tree used by the scheme
   * \param [in]     elements_index   The index in the target element array.
   * \param [in]     elements_from_index The index in the source element array.
   * \return                        The number of elements created in the target array.
   */
  template <adapt_action action>
  t8_locidx_t manipulate_elements (t8_element_array_t *elements,
                           const t8_element_array_t *const elements_from,
                           const t8_scheme *scheme,
                           const t8_eclass_t tree_class,
                           const t8_locidx_t elements_index,
                           const t8_locidx_t elements_from_index);

  template <>
  t8_locidx_t manipulate_elements<adapt_action::KEEP> (t8_element_array_t *elements,
                           const t8_element_array_t *const elements_from,
                           const t8_scheme *scheme,
                           const t8_eclass_t tree_class,
                           const t8_locidx_t elements_index,
                           const t8_locidx_t elements_from_index)
  {
    t8_element_t *element = t8_element_array_push (elements);
    const t8_element_t * element_from = t8_element_array_index_locidx (elements, elements_index);
    scheme->element_copy (tree_class, element_from, element);
    return 1;
  };

  template <>
  t8_locidx_t manipulate_elements<adapt_action::COARSEN> (t8_element_array_t *elements,
                           const t8_element_array_t *const elements_from,
                           const t8_scheme *scheme,
                           const t8_eclass_t tree_class,
                           const t8_locidx_t elements_index,
                           const t8_locidx_t elements_from_index)
  {
    t8_element_t *element = t8_element_array_push (elements);
    const t8_element_t * element_from = t8_element_array_index_locidx (elements, elements_index);
    T8_ASSERT (scheme->element_get_level (tree_class, element_from) > 0);
    scheme->element_get_parent (tree_class, element_from, element);

    /* Hier eventuell noch was mit num_children = num_siblings*/
    return 1;
  };

  template <>
  t8_locidx_t manipulate_elements<adapt_action::REFINE> (t8_element_array_t *elements,
                           const t8_element_array_t *const elements_from,
                           const t8_scheme *scheme,
                           const t8_eclass_t tree_class,
                           const t8_locidx_t elements_index,
                           const t8_locidx_t elements_from_index)
  {
    const t8_element_t * element_from = t8_element_array_index_locidx (elements_from, elements_from_index);
    const int num_children = scheme->element_get_num_children (tree_class, element_from);
    /* CONTINUE WORK HERE */
    (void) t8_element_array_push_count (elements, num_children);
    std::vector<t8_element_t *> children(num_children);
    for (int ichildren = 0; ichildren < num_children; ichildren++) {
      children[ichildren] = t8_element_array_index_locidx_mutable (elements, elements_index + ichildren);
    }
    scheme->element_get_children (tree_class, element_from, num_children, children.data());
    return num_children;
  };


  /**  * Class implementing a basic adaptation strategy for a forest of trees.
   */
class basic_adaptation {
  public:
    /** Constructor for basic_adaptation class.
     * \param [in] forest_in        The forest to be adapted.
     * \param [in] callback_in      The callback function to determine adaptation actions.
     */
    basic_adaptation (t8_forest_t forest, t8_forest_t forest_from, element_callback callback_in)
      : forest (forest), forest_from (forest_from), callback (callback_in) 
    {
      T8_ASSERT (forest != nullptr);
      T8_ASSERT (callback);
      if (forest_from != nullptr) {
        t8_forest_ref (forest_from);
      }
      T8_ASSERT (forest != nullptr);
      if (forest != nullptr) {
        t8_forest_ref (forest);
      }
    }

    /** Destructor for basic_adaptation class. */
    ~basic_adaptation () {
      if (forest_from != nullptr) {
        t8_forest_unref (&forest_from);
      }
      if (forest != nullptr) {
        t8_forest_unref (&forest);
      }
    }

    /** Perform the adaptation process on the forest. */
    void adapt();

    element_callback callback;                /**< The callback function to determine adaptation actions. */
  private:
    /**
     * Profile the adaptation process.
     */
    inline void 
    profile_adaptation(){
      T8_ASSERT(forest->profile != nullptr);
      forest->profile->adapt_runtime = -sc_MPI_Wtime();
    }

    /**
     * Collect adaptation actions for all elements in the source forest.
     */
    inline void
    collect_adapt_actions(){
      t8_locidx_t el_offset = 0;
      const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
      const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
      adapt_actions.resize (local_num_elements);

      const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

      /* For each element get the adaptation action */
      for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
          const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
          const t8_eclass_t tree_class = tree_from->eclass;
          const t8_element_array_t *elements_from = &tree_from->leaf_elements;
          const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (elements_from);
          for (t8_locidx_t i = 0; i < num_el_from; i++) {
              const t8_element_t * element = t8_element_array_index_locidx (elements_from, i);
              adapt_actions[el_offset + i] = callback (forest_from, ltree_id, element, scheme, tree_class);
          }
          el_offset += num_el_from;
      }
    };

    inline bool
    family_check(const t8_element_array_t *tree_elements_from, std::vector<const t8_element_t *> &elements_from, const t8_locidx_t offset, const t8_scheme *scheme, const t8_eclass_t tree_class){
      const int num_siblings = scheme->element_get_num_siblings (tree_class, elements_from[offset]);
      for (int isibling = 0; isibling < num_siblings; isibling++) {
        elements_from[isibling] = (const t8_element_t *) t8_element_array_index_locidx (tree_elements_from, offset + (t8_locidx_t )isibling);
        if (scheme->element_get_child_id (tree_class, elements_from[isibling]) != isibling) {
          return false;
        }
      }
      /* elements_are_family expects t8_element_t *const *; build a non-const pointer array */
      std::vector<t8_element_t *> children_nonconst(num_siblings);
      for (int i = 0; i < num_siblings; ++i)
        children_nonconst[i] = const_cast<t8_element_t *>(elements_from[i]);
      const bool is_family = scheme->elements_are_family (tree_class, children_nonconst.data());
      return is_family;
    }

    t8_forest_t forest;                       /**< The target forest */ 
    t8_forest_t forest_from;                  /**< The source forest to adapt from. */
    std::vector<adapt_action> adapt_actions;  /**< The adaptation actions for each element in the source forest. */
    bool profiling = false;                   /**< Flag to indicate if profiling is enabled. */
};

};

#endif /* T8_FOREST_ADAPT_HXX */
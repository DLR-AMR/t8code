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

/** \file t8_forest_adapt.hxx
 * CRTP base of the adaptation routine. The base drives the loop over all trees and elements of the source forest
 * and delegates every decision that differs between forest variants (complete forests, forests with holes,
 * subelements) to hooks, which a derived adaptor can hide to change the behavior.
 */

#ifndef T8_FOREST_ADAPT_HXX
#define T8_FOREST_ADAPT_HXX

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_data/t8_containers.h>

#include <array>

/** The actions an adapt callback \ref t8_forest_adapt_t can request for an element.
 * Return values greater than \ref t8_adapt_action::refine request a special refinement
 * (e.g. a subelement type) and are handled by \ref t8_adaptor::refine_special. */
enum struct t8_adapt_action : int {
  remove = -2,  /**< Remove the element. */
  coarsen = -1, /**< Coarsen the family into its parent. */
  keep = 0,     /**< Keep the element unchanged. */
  refine = 1    /**< Refine the element into its children. */
};

/** The state of the tree that is currently adapted. */
struct t8_adapt_tree_context
{
  const t8_locidx_t ltree_id;    /**< The local id of the tree. */
  const t8_eclass_t tree_class;  /**< The eclass of the tree. */
  const t8_scheme *const scheme; /**< The scheme of the forest. */
  /** The leaf elements of the tree in the source forest. They are not modified, the pointer is only mutable
   * since the adapt callback \ref t8_forest_adapt_t takes mutable elements. */
  t8_element_array_t *const elements_from;
  t8_element_array_t *const elements;  /**< The leaf elements of the tree in the new forest. */
  const t8_locidx_t num_elements_from; /**< The number of elements in \a elements_from. */
  t8_locidx_t el_considered;           /**< The index of the current element in \a elements_from. */
};

/** The elements passed to the adapt callback: either a family or a single element. */
struct t8_adapt_family
{
  /** The elements of the family. Only the first \a num_elements entries are defined. */
  std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN> elements;
  int num_elements; /**< The number of defined entries in \a elements. */
  bool is_family;   /**< True if \a elements form a family, false if \a elements is a single element. */
};

/** CRTP base of all adaptors. \ref adapt runs the fixed control flow and calls the hooks of \a TDerived.
 * The base implements every hook for complete forests. A derived adaptor hides a hook by declaring a member
 * function with the same name and signature, the base then calls the derived version without virtual dispatch.
 * \tparam TDerived The derived adaptor.
 */
template <typename TDerived>
struct t8_adaptor
{
  /** Constructor.
   * \param [in] forest_in The forest to be built by adapting its \a set_from forest. Its trees have to be
   *                       initialized without elements, see \ref t8_forest_copy_trees.
   */
  explicit t8_adaptor (t8_forest_t forest_in): forest (forest_in), forest_from (forest_in->set_from)
  {
    T8_ASSERT (forest != nullptr);
    T8_ASSERT (forest_from != nullptr);
    T8_ASSERT (forest->set_adapt_fn != nullptr);
  }

  /** Build the elements of the new forest by adapting the elements of the source forest. */
  void
  adapt ()
  {
    /* if profiling is enabled, measure runtime */
    if (forest->profile != nullptr) {
      forest->profile->adapt_runtime = -sc_MPI_Wtime ();
      /* DO NOT DELETE THE FOLLOWING line.
       * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
       * runtimes were computed to 0.
       * Only delete the line, if you know what you are doing. */
      t8_global_productionf ("Start adapt %f %f\n", sc_MPI_Wtime (), forest->profile->adapt_runtime);
    }
    t8_global_productionf ("Into t8_forest_adapt from %lld total elements\n",
                           (long long) forest_from->global_num_leaf_elements);

    T8_ASSERT (forest_from->incomplete_trees != -1);
    T8_ASSERT (forest->incomplete_trees == -1);
    T8_ASSERT (forest->trees->elem_count == forest_from->trees->elem_count);

    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);
    /* The number of elements in the new forest in all trees before the current one. */
    t8_locidx_t el_offset = 0;
    forest->local_num_leaf_elements = 0;
    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree = t8_forest_get_tree (forest, ltree_id);
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (&tree_from->leaf_elements);
      T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
      /* Elements cannot be inserted into an empty tree, so there is nothing to adapt. */
      if (num_el_from == 0) {
        continue;
      }
      T8_ASSERT (t8_element_array_get_count (&tree->leaf_elements) == 0);
      t8_adapt_tree_context context {
        ltree_id, tree->eclass, scheme, &tree_from->leaf_elements, &tree->leaf_elements, num_el_from, 0
      };
      while (context.el_considered < num_el_from) {
        adapt_family (context);
      }

      const t8_locidx_t el_inserted = (t8_locidx_t) t8_element_array_get_count (&tree->leaf_elements);
      /* It is not supported to delete all elements from a tree. In this case, we abort. */
      SC_CHECK_ABORTF (el_inserted != 0,
                       "ERROR: All elements of tree %i were removed. Removing all elements of a tree "
                       "is currently not supported. See also https://github.com/DLR-AMR/t8code/issues/1137.",
                       ltree_id);
      tree->elements_offset = el_offset;
      el_offset += el_inserted;
      forest->local_num_leaf_elements += el_inserted;
      /* Release the memory that is not needed anymore. */
      t8_element_array_resize (&tree->leaf_elements, el_inserted);
    }

    t8_forest_comm_global_num_leaf_elements (forest);
    derived ().finalize ();

    t8_global_productionf ("Done t8_forest_adapt with %lld total elements\n",
                           (long long) forest->global_num_leaf_elements);
    /* if profiling is enabled, measure runtime */
    if (forest->profile != nullptr) {
      forest->profile->adapt_runtime += sc_MPI_Wtime ();
      /* DO NOT DELETE THE FOLLOWING line.
       * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
       * runtimes were computed to 0.
       * Only delete the line, if you know what you are doing. */
      t8_global_productionf ("End adapt %f %f\n", sc_MPI_Wtime (), forest->profile->adapt_runtime);
    }
  }

  /** Collect the element at \a context.el_considered and, if they form a family, its siblings.
   * \param [in] context The current tree.
   * \return The family starting at \a context.el_considered, or this element alone.
   */
  t8_adapt_family
  get_family (const t8_adapt_tree_context &context)
  {
    t8_adapt_family family;
    const t8_element_t *first = t8_element_array_index_locidx (context.elements_from, context.el_considered);
    const int num_siblings = context.scheme->element_get_num_siblings (context.tree_class, first);
    T8_ASSERT (num_siblings <= T8_ECLASS_MAX_CHILDREN);
    /* The number of consecutive elements whose child ids are 0, 1, ..., num_collected - 1. */
    int num_collected = 0;
    for (; num_collected < num_siblings && context.el_considered + num_collected < context.num_elements_from;
         num_collected++) {
      family.elements[num_collected]
        = t8_element_array_index_locidx_mutable (context.elements_from, context.el_considered + num_collected);
      /* Quick check: a family has the child ids 0, 1, ..., num_siblings - 1. This is not sufficient, therefore
       * we check with elements_are_family below. */
      if (context.scheme->element_get_child_id (context.tree_class, family.elements[num_collected]) != num_collected) {
        break;
      }
    }
    family.is_family = num_collected == num_siblings
                       && context.scheme->elements_are_family (context.tree_class, family.elements.data ());
    family.num_elements = family.is_family ? num_siblings : 1;
    return family;
  }

  /** Ask the adapt callback of the forest what to do with \a family.
   * \param [in] context The current tree.
   * \param [in] family  The elements to pass to the callback.
   * \return The return value of the callback, see \ref t8_forest_adapt_t.
   */
  int
  query (const t8_adapt_tree_context &context, t8_adapt_family &family)
  {
    return forest->set_adapt_fn (forest, forest_from, context.ltree_id, context.tree_class, context.el_considered,
                                 context.scheme, family.is_family, family.num_elements, family.elements.data ());
  }

  /** Insert the children of the first element of \a family into the new tree.
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  refine (t8_adapt_tree_context &context, const t8_adapt_family &family)
  {
    const int num_children = context.scheme->element_get_num_children (context.tree_class, family.elements[0]);
    (void) t8_element_array_push_count (context.elements, num_children);
    context.scheme->element_get_children (context.tree_class, family.elements[0], num_children,
                                          last_elements (context, num_children).data ());
    return 1;
  }

  /** Replace \a family by its parent in the new tree.
   * \param [in, out] context The current tree.
   * \param [in]      family  The family to coarsen.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  coarsen (t8_adapt_tree_context &context, const t8_adapt_family &family)
  {
    T8_ASSERT (family.is_family);
    T8_ASSERT (context.scheme->element_get_level (context.tree_class, family.elements[0]) > 0);
    t8_element_t *parent = t8_element_array_push (context.elements);
    context.scheme->element_get_parent (context.tree_class, family.elements[0], parent);
    return family.num_elements;
  }

  /** Copy the first element of \a family into the new tree.
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  keep (t8_adapt_tree_context &context, const t8_adapt_family &family)
  {
    t8_element_t *element = t8_element_array_push (context.elements);
    context.scheme->element_copy (context.tree_class, family.elements[0], element);
    return 1;
  }

  /** Drop the first element of \a family from the new tree.
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  remove ([[maybe_unused]] t8_adapt_tree_context &context, [[maybe_unused]] const t8_adapt_family &family)
  {
    element_removed = true;
    return 1;
  }

  /** Handle a callback return value greater than \ref t8_adapt_action::refine.
   * The base does not support any special refinement.
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \param [in]      type    The return value of the callback minus one.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  refine_special ([[maybe_unused]] t8_adapt_tree_context &context, [[maybe_unused]] const t8_adapt_family &family,
                  const int type)
  {
    SC_ABORTF ("The adapt callback returned %i, which is not supported by this forest.", type + 1);
    return 0;
  }

  /** Called after all trees are adapted and the global number of elements is known.
   * Determines whether the new forest has incomplete trees.
   */
  void
  finalize ()
  {
    /* If the old forest already contained incomplete trees, the new one does too. */
    if (forest_from->incomplete_trees) {
      T8_ASSERT (forest_from->incomplete_trees == 1);
      forest->incomplete_trees = 1;
      return;
    }
    /* Inform the other processes whether elements were removed locally. */
    const int local_removed = element_removed;
    int incomplete_trees;
    const int mpiret = sc_MPI_Allreduce (&local_removed, &incomplete_trees, 1, sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (incomplete_trees == 1 || incomplete_trees == 0);
    forest->incomplete_trees = incomplete_trees;
  }

 protected:
  /** Get pointers to the last \a count elements of the new tree.
   * \param [in] context The current tree.
   * \param [in] count   The number of elements.
   * \return The pointers to the last \a count elements, the other entries are undefined.
   */
  static std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN>
  last_elements (const t8_adapt_tree_context &context, const int count)
  {
    T8_ASSERT (count <= T8_ECLASS_MAX_CHILDREN);
    std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN> last;
    const t8_locidx_t first = (t8_locidx_t) t8_element_array_get_count (context.elements) - count;
    for (int ielement = 0; ielement < count; ielement++) {
      last[ielement] = t8_element_array_index_locidx_mutable (context.elements, first + ielement);
    }
    return last;
  }

  t8_forest_t forest;           /**< The forest under construction. */
  t8_forest_t forest_from;      /**< The forest that is adapted. */
  bool element_removed = false; /**< True if an element was removed on this process, false otherwise. */

 private:
  /** Query the callback for the element at \a context.el_considered and apply the requested action.
   * \param [in, out] context The current tree, \a context.el_considered is advanced past the consumed elements.
   */
  void
  adapt_family (t8_adapt_tree_context &context)
  {
    t8_adapt_family family = derived ().get_family (context);
    T8_ASSERT (family.num_elements > 0);
    int refine_value = derived ().query (context, family);
    T8_ASSERT (family.is_family || refine_value != static_cast<int> (t8_adapt_action::coarsen));
    if (refine_value > 0
        && (context.scheme->element_get_level (context.tree_class, family.elements[0]) >= forest->maxlevel
            || !context.scheme->element_is_refinable (context.tree_class, family.elements[0]))) {
      /* Only refine an element if it does not exceed the maximum level and if it is refinable */
      refine_value = static_cast<int> (t8_adapt_action::keep);
    }

    if (refine_value > static_cast<int> (t8_adapt_action::refine)) {
      context.el_considered += derived ().refine_special (context, family, refine_value - 1);
      return;
    }
    switch (static_cast<t8_adapt_action> (refine_value)) {
    case t8_adapt_action::refine:
      context.el_considered += derived ().refine (context, family);
      break;
    case t8_adapt_action::coarsen:
      context.el_considered += derived ().coarsen (context, family);
      break;
    case t8_adapt_action::keep:
      context.el_considered += derived ().keep (context, family);
      break;
    case t8_adapt_action::remove:
      context.el_considered += derived ().remove (context, family);
      break;
    default:
      SC_ABORTF ("The adapt callback returned the invalid value %i.", refine_value);
    }
  }

  /** \return This object as the derived adaptor. */
  TDerived &
  derived ()
  {
    return static_cast<TDerived &> (*this);
  }
};

#endif /* T8_FOREST_ADAPT_HXX */

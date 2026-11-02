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

#include <algorithm>
#include <array>
#include <vector>

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
  /** The index of the first element in \a elements that may be coarsened recursively.
   * Elements before it stem from a recursive refinement and must not be coarsened again. */
  t8_locidx_t el_coarsen;
};

/** The elements passed to the adapt callback: either a family or a single element. */
struct t8_adapt_family
{
  /** The elements of the family. Only the first \a num_elements entries are defined. */
  std::array<t8_element_t *, T8_ECLASS_MAX_CHILDREN> elements;
  int num_elements; /**< The number of defined entries in \a elements. */
  bool is_family;   /**< True if \a elements form a family, false if \a elements is a single element. */
};

/** Callback policy of \ref t8_adaptor that asks the per-element adapt callback \ref t8_forest_adapt_t of the forest
 * for each element or family.
 */
struct t8_adapt_element_callback
{
  /** Constructor.
   * \param [in] forest The forest under construction. Its adapt callback has to be set.
   */
  explicit t8_adapt_element_callback ([[maybe_unused]] const t8_forest_t forest)
  {
    T8_ASSERT (forest->set_adapt_fn != nullptr);
  }

  /** Called before the elements of a tree are adapted. The per-element callback needs no preparation.
   * \param [in] forest      The forest under construction.
   * \param [in] forest_from The forest that is adapted.
   * \param [in] context     The tree that is adapted next.
   */
  void
  begin_tree ([[maybe_unused]] const t8_forest_t forest, [[maybe_unused]] const t8_forest_t forest_from,
              [[maybe_unused]] const t8_adapt_tree_context &context)
  {
  }

  /** Ask the adapt callback of the forest what to do with \a family.
   * \param [in] forest      The forest under construction.
   * \param [in] forest_from The forest that is adapted.
   * \param [in] context     The current tree.
   * \param [in] family      The elements to pass to the callback.
   * \return The return value of the callback, see \ref t8_forest_adapt_t.
   */
  int
  query (const t8_forest_t forest, const t8_forest_t forest_from, const t8_adapt_tree_context &context,
         t8_adapt_family &family)
  {
    return forest->set_adapt_fn (forest, forest_from, context.ltree_id, context.tree_class, context.el_considered,
                                 context.scheme, family.is_family, family.num_elements, family.elements.data ());
  }
};

/** Callback policy of \ref t8_adaptor that asks the batched adapt callback \ref t8_forest_adapt_batched_t of the
 * forest once per tree and answers the queries for the elements and families of the tree from its result.
 * \note The actions refer to the elements of the source forest, so this policy cannot be used for recursive
 *       adaptation.
 */
struct t8_adapt_batched_callback
{
  /** Constructor.
   * \param [in] forest The forest under construction. Its batched adapt callback has to be set.
   */
  explicit t8_adapt_batched_callback ([[maybe_unused]] const t8_forest_t forest)
  {
    T8_ASSERT (forest->set_adapt_batched_fn != nullptr);
  }

  /** Call the batched adapt callback for all elements of the tree.
   * \param [in] forest      The forest under construction.
   * \param [in] forest_from The forest that is adapted.
   * \param [in] context     The tree that is adapted next.
   */
  void
  begin_tree (const t8_forest_t forest, const t8_forest_t forest_from, const t8_adapt_tree_context &context)
  {
    actions.assign ((size_t) context.num_elements_from, static_cast<int> (t8_adapt_action::keep));
    forest->set_adapt_batched_fn (forest, forest_from, context.ltree_id, context.tree_class, context.scheme,
                                  context.elements_from, actions.data ());
  }

  /** Get the action for \a family from the actions of the tree. A family is coarsened if all its members are
   * marked for coarsening. Otherwise, a coarsening mark of the first element is treated as keep.
   * \param [in] forest      The forest under construction.
   * \param [in] forest_from The forest that is adapted.
   * \param [in] context     The current tree.
   * \param [in] family      The elements starting at \a context.el_considered.
   * \return The adapt action, see \ref t8_forest_adapt_t.
   */
  int
  query ([[maybe_unused]] const t8_forest_t forest, [[maybe_unused]] const t8_forest_t forest_from,
         const t8_adapt_tree_context &context, const t8_adapt_family &family) const
  {
    T8_ASSERT (context.el_considered + family.num_elements <= (t8_locidx_t) actions.size ());
    const int coarsen = static_cast<int> (t8_adapt_action::coarsen);
    const auto first = actions.begin () + context.el_considered;
    if (family.is_family && std::all_of (first, first + family.num_elements, [coarsen] (const int action) {
          return action == coarsen;
        })) {
      return coarsen;
    }
    return *first == coarsen ? static_cast<int> (t8_adapt_action::keep) : *first;
  }

  std::vector<int> actions; /**< The actions of the elements of the current tree, as set by the callback. */
};

/** CRTP base of all adaptors. \ref adapt runs the fixed control flow and calls the hooks of \a TDerived.
 * The base implements every hook for complete forests. A derived adaptor hides a hook by declaring a member
 * function with the same name and signature, the base then calls the derived version without virtual dispatch.
 * \tparam TDerived   The derived adaptor.
 * \tparam TRecursive If true, newly created elements are passed to the callback again and refined or coarsened
 *                    recursively, see \ref t8_forest_set_adapt.
 * \tparam TCallback  The policy that asks the user callback for the adapt actions, e.g.
 *                    \ref t8_adapt_element_callback.
 */
template <typename TDerived, bool TRecursive, typename TCallback>
struct t8_adaptor
{
  /** Constructor.
   * \param [in] forest_in The forest to be built by adapting its \a set_from forest. Its trees have to be
   *                       initialized without elements, see \ref t8_forest_copy_trees.
   */
  explicit t8_adaptor (t8_forest_t forest_in)
    : forest (forest_in), forest_from (forest_in->set_from), callback (forest_in)
  {
    T8_ASSERT (forest != nullptr);
    T8_ASSERT (forest_from != nullptr);
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
        ltree_id, tree->eclass, scheme, &tree_from->leaf_elements, &tree->leaf_elements, num_el_from, 0, 0
      };
      callback.begin_tree (forest, forest_from, context);
      while (context.el_considered < num_el_from) {
        adapt_family (context);
      }
      if constexpr (TRecursive) {
        destroy_refine_buffer (context);
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
    return collect_family (context, context.elements_from, context.el_considered, context.num_elements_from);
  }

  /** Collect the elements of the new tree from \a first to its end for recursive coarsening.
   * \param [in] context The current tree.
   * \param [in] first   The index of the first possible family member in \a context.elements.
   * \return The family starting at \a first, or this element alone.
   */
  t8_adapt_family
  get_new_family (const t8_adapt_tree_context &context, const t8_locidx_t first)
  {
    return collect_family (context, context.elements, first,
                           (t8_locidx_t) t8_element_array_get_count (context.elements));
  }

  /** Ask the callback policy what to do with \a family.
   * \param [in] context The current tree.
   * \param [in] family  The elements to decide on.
   * \return The adapt action, see \ref t8_forest_adapt_t.
   */
  int
  query (const t8_adapt_tree_context &context, t8_adapt_family &family)
  {
    return callback.query (forest, forest_from, context, family);
  }

  /** Insert the children of the first element of \a family into the new tree.
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  refine (t8_adapt_tree_context &context, const t8_adapt_family &family)
  {
    if constexpr (TRecursive) {
      refine_recursive (context, family.elements[0]);
      context.el_coarsen = (t8_locidx_t) t8_element_array_get_count (context.elements);
      return 1;
    }
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
    if constexpr (TRecursive) {
      coarsen_recursive (context);
    }
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
    if constexpr (TRecursive) {
      coarsen_recursive (context);
    }
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

  /** Collect the element at \a first in \a array and, if they form a family, its siblings.
   * \param [in] context The current tree.
   * \param [in] array   Either \a context.elements_from or \a context.elements.
   * \param [in] first   The index of the first possible family member in \a array.
   * \param [in] end     The number of elements in \a array that may be part of the family.
   * \return The family starting at \a first, or this element alone.
   */
  static t8_adapt_family
  collect_family (const t8_adapt_tree_context &context, t8_element_array_t *array, const t8_locidx_t first,
                  const t8_locidx_t end)
  {
    t8_adapt_family family;
    const int num_siblings
      = context.scheme->element_get_num_siblings (context.tree_class, t8_element_array_index_locidx (array, first));
    T8_ASSERT (num_siblings <= T8_ECLASS_MAX_CHILDREN);
    /* The number of consecutive elements whose child ids are 0, 1, ..., num_collected - 1. */
    int num_collected = 0;
    for (; num_collected < num_siblings && first + num_collected < end; num_collected++) {
      family.elements[num_collected] = t8_element_array_index_locidx_mutable (array, first + num_collected);
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

  /** Refine \a parent and pass each child to the callback, which may refine it further.
   * The resulting leaves are appended to the new tree in SFC order.
   * \param [in, out] context The current tree.
   * \param [in]      parent  The element to refine. Its level is smaller than the maximum level of the forest.
   */
  void
  refine_recursive (t8_adapt_tree_context &context, const t8_element_t *parent)
    requires (TRecursive)
  {
    const int child_level = context.scheme->element_get_level (context.tree_class, parent) + 1;
    T8_ASSERT (child_level <= forest->maxlevel);
    const int max_children = context.scheme->get_max_num_children (context.tree_class);
    if (refine_buffer.empty ()) {
      refine_buffer.resize ((size_t) (forest->maxlevel + 1) * max_children);
      context.scheme->element_new (context.tree_class, (int) refine_buffer.size (), refine_buffer.data ());
    }
    /* The children of all elements of one level share a block of the buffer, since we refine depth first. */
    t8_element_t **children = refine_buffer.data () + (size_t) child_level * max_children;
    const int num_children = context.scheme->element_get_num_children (context.tree_class, parent);
    context.scheme->element_get_children (context.tree_class, parent, num_children, children);

    for (int ichild = 0; ichild < num_children; ichild++) {
      t8_adapt_family child { { children[ichild] }, 1, false };
      const int refine_value = derived ().query (context, child);
      T8_ASSERT (refine_value != static_cast<int> (t8_adapt_action::coarsen));
      if (refine_value == static_cast<int> (t8_adapt_action::refine) && child_level < forest->maxlevel
          && context.scheme->element_is_refinable (context.tree_class, children[ichild])) {
        refine_recursive (context, children[ichild]);
      }
      else if (refine_value == static_cast<int> (t8_adapt_action::remove)) {
        element_removed = true;
      }
      else {
        t8_element_t *element = t8_element_array_push (context.elements);
        context.scheme->element_copy (context.tree_class, children[ichild], element);
      }
    }
  }

  /** Free the buffer used by \ref refine_recursive, since the next tree may have a different eclass.
   * \param [in] context The tree that was adapted.
   */
  void
  destroy_refine_buffer (const t8_adapt_tree_context &context)
    requires (TRecursive)
  {
    if (!refine_buffer.empty ()) {
      context.scheme->element_destroy (context.tree_class, (int) refine_buffer.size (), refine_buffer.data ());
      refine_buffer.clear ();
    }
  }

  /** If the last element of the new tree is the last member of a family, pass the family to the callback and
   * coarsen it if requested. Repeat with the resulting parent.
   * \param [in, out] context The current tree.
   */
  void
  coarsen_recursive (t8_adapt_tree_context &context)
    requires (TRecursive)
  {
    t8_locidx_t first = find_family_start (context, (t8_locidx_t) t8_element_array_get_count (context.elements) - 1);
    while (first >= context.el_coarsen) {
      t8_adapt_family family = derived ().get_new_family (context, first);
      if (!family.is_family || derived ().query (context, family) != static_cast<int> (t8_adapt_action::coarsen)) {
        return;
      }
      T8_ASSERT (context.scheme->element_get_level (context.tree_class, family.elements[0]) > 0);
      context.scheme->element_get_parent (context.tree_class, family.elements[0], family.elements[0]);
      /* The parent replaces the family. */
      t8_element_array_resize (context.elements, first + 1);
      first = find_family_start (context, first);
    }
  }

  /** Find the first member of the family whose last member is at \a last in the new tree.
   * \param [in] context The current tree.
   * \param [in] last    The index of an element in \a context.elements.
   * \return The index of the first family member, or a negative value if the element at \a last is not the
   *         last member of a family that can be coarsened.
   */
  t8_locidx_t
  find_family_start (const t8_adapt_tree_context &context, const t8_locidx_t last)
    requires (TRecursive)
  {
    const t8_scheme *scheme = context.scheme;
    const t8_eclass_t tree_class = context.tree_class;
    const t8_locidx_t num_elements = (t8_locidx_t) t8_element_array_get_count (context.elements);
    T8_ASSERT (0 <= last && last < num_elements);
    const t8_element_t *element = t8_element_array_index_locidx (context.elements, last);
    const int level = scheme->element_get_level (tree_class, element);
    const int num_siblings = scheme->element_get_num_siblings (tree_class, element);
    const int child_id = scheme->element_get_child_id (tree_class, element);
    /* Elements whose family consists of only one element (vertices) and elements of level 0 cannot be coarsened.
     * Otherwise, the element has to be the last member of its family. */
    if (!(child_id > 0 && child_id == num_siblings - 1) || level == 0) {
      return -1;
    }

    t8_element_t *parent;
    t8_element_t *parent_compare;
    scheme->element_new (tree_class, 1, &parent);
    scheme->element_new (tree_class, 1, &parent_compare);
    scheme->element_get_parent (tree_class, element, parent);

    /* Loop backward over all possible family members until we hit an element that is not part of the family
     * or we have reached the maximum number of members. */
    t8_locidx_t iback = 1;
    t8_locidx_t first = last;
    const t8_element_t *element_compare = nullptr;
    for (; iback < num_siblings && iback < num_elements; iback++) {
      element_compare = t8_element_array_index_locidx (context.elements, last - iback);
      /* Comparing the levels first may avoid the more expensive comparison of the parents. */
      if (scheme->element_get_level (tree_class, element_compare) != level) {
        break;
      }
      scheme->element_get_parent (tree_class, element_compare, parent_compare);
      if (!scheme->element_is_equal (tree_class, parent, parent_compare)) {
        break;
      }
      first = last - iback;
    }

    /* If we found fewer elements than a family has, the element before them may be a descendant of the parent.
     * Then coarsening would overlap it. */
    if (iback < num_siblings && iback < num_elements
        && scheme->element_get_level (tree_class, element_compare) > level) {
      scheme->element_get_nca (tree_class, element, element_compare, parent_compare);
      T8_ASSERT (scheme->element_get_level (tree_class, parent_compare) <= level - 1);
      if (scheme->element_get_level (tree_class, parent_compare) == level - 1) {
        first = -1;
      }
    }
    scheme->element_destroy (tree_class, 1, &parent);
    scheme->element_destroy (tree_class, 1, &parent_compare);

#if T8_ENABLE_MPI
    /* The first element on process rank must have child_id 0, otherwise other family members could be on
     * process rank - 1. */
    if (first == 0 && forest->mpirank > 0
        && scheme->element_get_child_id (tree_class, t8_element_array_index_locidx (context.elements, 0)) > 0) {
      return -1;
    }
#endif
    return first;
  }

  t8_forest_t forest;           /**< The forest under construction. */
  t8_forest_t forest_from;      /**< The forest that is adapted. */
  TCallback callback;           /**< Asks the user callback for the adapt actions. */
  bool element_removed = false; /**< True if an element was removed on this process, false otherwise. */
  /** The elements used by \ref refine_recursive, one block of the maximum number of children per level.
   * Empty if no element of the current tree was refined recursively. */
  std::vector<t8_element_t *> refine_buffer;

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

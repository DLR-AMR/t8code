/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_forest_iterate.cxx
 * Implements functions declared in \ref t8_forest_iterate.h.
 */

#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>

#include <algorithm>

void
t8_forest_split_array (const t8_element_t *element, const t8_element_array_view &leaf_elements, size_t *offsets)
{
  const t8_scheme *scheme = leaf_elements.get_scheme ();
  const t8_eclass_t tree_class = leaf_elements.get_tree_class ();
  const int num_children = scheme->element_get_num_children (tree_class, element);
  /** The refinement level of the children of \a element. */
  const int child_level = scheme->element_get_level (tree_class, element) + 1;

  /* The leaves are sorted, so the descendants of each child of element are stored contiguously and ordered by the
   * child id, which is the leaves' ancestor id at the level of the children. In other words for each child C of
   * element, find the indices i, j such that all descendants of C are elements[i], ..., elements[j-1].
   * We binary search the first leaf of each child. Every probed leaf with child id c is an upper bound for the
   * first leaf of all children up to c, so we store it in their offsets. This narrows the search range of the
   * following children, as in sc_array_split. Until its search is finished, offsets[ichild] is an upper bound. */
  const size_t num_leaves = leaf_elements.size ();
  offsets[0] = 0;
  std::fill (offsets + 1, offsets + num_children + 1, num_leaves);
  /** All leaves before this index are descendants of children with a smaller id than the currently searched child. */
  size_t low = 0;
  for (int ichild = 1; ichild < num_children; ++ichild) {
    while (low < offsets[ichild]) {
      const size_t guess = low + (offsets[ichild] - low) / 2;
      T8_ASSERT (child_level <= scheme->element_get_level (tree_class, leaf_elements[guess]));
      const int guess_child_id = scheme->element_get_ancestor_id (tree_class, leaf_elements[guess], child_level);
      T8_ASSERT (0 <= guess_child_id && guess_child_id < num_children);
      if (guess_child_id < ichild) {
        low = guess + 1;
      }
      else {
        std::fill (offsets + ichild, offsets + guess_child_id + 1, guess);
      }
    }
  }
}

void
t8_forest_iterate_faces (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                         const int face, const t8_element_array_view &leaf_elements,
                         const t8_locidx_t tree_lindex_of_first_leaf, const t8_forest_iterate_face_fn callback,
                         void *user_data)
{
  t8_debugf ("Entering t8_forest_iterate_faces with leaf_index %i and %zu total leaves.\n", tree_lindex_of_first_leaf,
             leaf_elements.size ());
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
#if T8_ENABLE_DEBUG
  const t8_locidx_t num_ghost_trees = t8_forest_get_num_ghost_trees (forest);
  const t8_locidx_t num_local_and_ghost_trees = num_local_trees + num_ghost_trees;
  T8_ASSERT (0 <= ltreeid && ltreeid < num_local_and_ghost_trees);
#endif

  // Check whether we are in a local tree or ghost tree
  const bool tree_is_local = t8_forest_tree_is_local (forest, ltreeid);
  const bool tree_is_ghost = !tree_is_local;
  const t8_locidx_t local_or_ghost_tree_id = tree_is_local ? ltreeid : ltreeid - num_local_trees;

  const size_t elem_count = leaf_elements.size ();
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }
  const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  // TODO: This does a costly binary search. Is there a better criterion to check
  //        whether element is a leaf or not?
  //        We tried (elem_count == 1) but this does not work since we could
  //        start the search with a full family and element being one sibling.
  //        In that case, element is a leaf but elem_count is not 1.
  bool is_leaf = t8_forest_element_is_leaf_or_ghost (forest, element, local_or_ghost_tree_id, tree_is_ghost);

#if T8_ENABLE_DEBUG
  if (!is_leaf) {
    /* Check whether element has smaller level than the first leaf */
    const t8_element_t *leaf = leaf_elements[0];
    T8_ASSERT (t8_forest_element_is_leaf_or_ghost (forest, leaf, local_or_ghost_tree_id, tree_is_ghost));
    T8_ASSERT (scheme->element_get_level (eclass, element) < scheme->element_get_level (eclass, leaf));

    // Verify that all leaves in leaf_elements are descendants of element
    for (const t8_element_t *ileaf : leaf_elements) {
      T8_ASSERT (scheme->element_is_ancestor (eclass, element, ileaf));
    }
  }
#endif

  /* Call the callback function element, if it returns true, we continue with the top-down recursion */
  const int ret
    = callback (forest, ltreeid, element, face, is_leaf, leaf_elements, tree_lindex_of_first_leaf, user_data);
  if (!ret || is_leaf) {
    // The callback returned false or the element is a leaf.
    // We abort the recursion.
    return;
  }

  /* Enter the recursion */
  /* We compute all face children of E, compute their leaf arrays and call iterate_faces */
  /* allocate the memory to store the face children */
  const int num_face_children = scheme->element_get_num_face_children (eclass, element, face);
  t8_element_t **face_children = T8_ALLOC (t8_element_t *, num_face_children);
  scheme->element_new (eclass, num_face_children, face_children);
  /* Memory for the child indices of the face children */
  int *child_indices = T8_ALLOC (int, num_face_children);
  /* Memory for the indices that split the leaf_elements array */
  const int num_children = scheme->element_get_num_children (eclass, element);
  size_t *split_offsets = T8_ALLOC (size_t, num_children + 1);
  /* Compute the face children */
  scheme->element_get_children_at_face (eclass, element, face, face_children, num_face_children, child_indices);
  /* Split the leaves array in portions belonging to the children of element */
  t8_forest_split_array (element, leaf_elements, split_offsets);
  for (int iface = 0; iface < num_face_children; iface++) {
    /* Check if there are any leaf elements for this face child */
    T8_ASSERT (0 <= child_indices[iface]);
    T8_ASSERT (child_indices[iface] < num_children + 1);
    const size_t indexa = split_offsets[child_indices[iface]];     /* first leaf of this face child */
    const size_t indexb = split_offsets[child_indices[iface] + 1]; /* first leaf of next child */
    t8_debugf ("Computed indices for face child %i:  %zd %zd\n", iface, indexa, indexb);
    if (indexa < indexb) {
      /* Compute the corresponding face number of this face child */
      const int child_face = scheme->element_face_get_child_face (eclass, element, face, iface);
      /* There exist leaves of this face child in leaf_elements,
       * we enter the recursion with a view of these leaves */
      t8_forest_iterate_faces (forest, ltreeid, face_children[iface], child_face,
                               leaf_elements.subview (indexa, indexb - indexa), indexa + tree_lindex_of_first_leaf,
                               callback, user_data);
    }
  }
  /* clean-up */
  scheme->element_destroy (eclass, num_face_children, face_children);
  T8_FREE (face_children);
  T8_FREE (child_indices);
  T8_FREE (split_offsets);
}

/* We want to export the remaining implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

void
t8_forest_iterate_replace (t8_forest_t forest_new, t8_forest_t forest_old, t8_forest_replace_t replace_fn)
{
  t8_global_productionf ("Into t8_forest_iterate_replace\n");
  T8_ASSERT (t8_forest_is_committed (forest_old));
  T8_ASSERT (t8_forest_is_committed (forest_new));
  const t8_scheme *scheme = t8_forest_get_scheme (forest_new);
  // Check that the two forests use the same scheme.
  T8_ASSERT (scheme == t8_forest_get_scheme (forest_old));

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest_new);
  T8_ASSERT (num_local_trees == t8_forest_get_num_local_trees (forest_old));

  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    /* Loop over the trees */
    /* Get the number of elements of this tree in old and new forest */
    const t8_locidx_t elems_per_tree_new = t8_forest_get_tree_num_leaf_elements (forest_new, itree);
    const t8_locidx_t elems_per_tree_old = t8_forest_get_tree_num_leaf_elements (forest_old, itree);
    /* Get the eclass of the tree */
    t8_eclass_t const tree_class = t8_forest_get_tree_class (forest_new, itree);
    T8_ASSERT (tree_class == t8_forest_get_tree_class (forest_old, itree));
    t8_locidx_t ielem_new = 0;
    t8_locidx_t ielem_old = 0;
    while (ielem_new < elems_per_tree_new) {
      /* Iterate over the elements */
      T8_ASSERT (ielem_new < elems_per_tree_new);
      T8_ASSERT (ielem_old < elems_per_tree_old);

      /* Get pointers to the elements */
      const t8_element_t *elem_new = t8_forest_get_leaf_element_in_tree (forest_new, itree, ielem_new);
      const t8_element_t *elem_old = t8_forest_get_leaf_element_in_tree (forest_old, itree, ielem_old);

      /* Get the levels of these elements */
      const int level_new = scheme->element_get_level (tree_class, elem_new);
      const int level_old = scheme->element_get_level (tree_class, elem_old);

      if (forest_new->incomplete_trees) {
        /* If el_removed is 1, the element in forest_new has been removed.
         * It is assumed that no element was removed. */
        int el_removed = 0;
        if (level_old < level_new) {
          /* elem_old got refined or removed */
          t8_element_t *elem_parent;
          scheme->element_new (tree_class, 1, &elem_parent);
          scheme->element_get_parent (tree_class, elem_new, elem_parent);
          if (scheme->element_is_equal (tree_class, elem_old, elem_parent)) {
            /* elem_old got refined */
            T8_ASSERT (level_new == level_old + 1);
            const t8_locidx_t family_size = scheme->element_get_num_children (tree_class, elem_old);
#if T8_DEBUG
            /* Check if family of new refined elements is complete */
            T8_ASSERT (ielem_new + family_size <= elems_per_tree_new);
            for (t8_locidx_t ielem = 1; ielem < family_size; ielem++) {
              const t8_element_t *elem_new_debug
                = t8_forest_get_leaf_element_in_tree (forest_new, itree, ielem_new + ielem);
              scheme->t8_element_parent (elem_new_debug, elem_parent);
              SC_CHECK_ABORT (scheme->element_is_equal (tree_class, elem_old, elem_parent), "Family is not complete.");
            }
#endif
            scheme->element_destroy (tree_class, 1, &elem_parent);
            const int refine = 1;
            replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, family_size,
                        ielem_new);
            /* Advance to the next element */
            ielem_new += family_size;
            ielem_old++;
          }
          else {
            /* elem_old got removed */
            el_removed = 1;
            scheme->element_destroy (tree_class, 1, &elem_parent);
          }
        }
        else if (level_old > level_new) {
          /* elem_old got coarsened or removed */
          t8_element_t *elem_parent;
          scheme->element_new (tree_class, 1, &elem_parent);
          scheme->element_get_parent (tree_class, elem_old, elem_parent);
          if (scheme->element_is_equal (tree_class, elem_new, elem_parent)) {
            /* elem_old got coarsened */
            T8_ASSERT (level_new == level_old - 1);
            /* Get size of family of old forest */
            int family_size = 1;
            for (t8_locidx_t ielem = 1; ielem < scheme->element_get_num_children (tree_class, elem_new)
                                        && ielem + ielem_old < elems_per_tree_old;
                 ielem++) {
              elem_old = t8_forest_get_leaf_element_in_tree (forest_old, itree, ielem_old + ielem);
              scheme->element_get_parent (tree_class, elem_old, elem_parent);
              if (scheme->element_is_equal (tree_class, elem_new, elem_parent)) {
                family_size++;
              }
            }
            T8_ASSERT (family_size <= scheme->element_get_num_children (tree_class, elem_new));
#if T8_ENABLE_DEBUG
            /* Check whether elem_old is the first element of the family */
            for (t8_locidx_t ielem = 1;
                 ielem < scheme->element_get_num_children (tree_class, elem_old) && ielem_old - ielem >= 0; ielem++) {
              const t8_element_t *elem_old_debug
                = t8_forest_get_leaf_element_in_tree (forest_old, itree, ielem_old - ielem);
              scheme->element_get_parent (tree_class, elem_old_debug, elem_parent);
              SC_CHECK_ABORT (!scheme->element_is_equal (tree_class, elem_new, elem_parent),
                              "elem_old is not the first of the family.");
            }
#endif
            scheme->element_destroy (tree_class, 1, &elem_parent);
            const int refine = -1;
            replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, family_size, ielem_old, 1,
                        ielem_new);
            /* Advance to the next element */
            ielem_new++;
            ielem_old += family_size;
          }
          else {
            /* elem_old got removed */
            el_removed = 1;
            scheme->element_destroy (tree_class, 1, &elem_parent);
          }
        }
        else {
          /* elem_old was untouched or got removed */
          if (scheme->element_is_equal (tree_class, elem_new, elem_old)) {
            /* elem_new = elem_old */
            const int refine = 0;
            replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, 1, ielem_new);
            /* Advance to the next element */
            ielem_new++;
            ielem_old++;
          }
          else {
            /* elem_old got removed */
            el_removed = 1;
          }
        }
        if (el_removed) {
          T8_ASSERT (el_removed == 1);
          T8_ASSERT (forest_new->incomplete_trees == 1);
          /* element got removed */
          const int refine = -2;
          replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, 0, -1);
          /* Advance to the next element */
          ielem_old++;
        }
        T8_ASSERT (el_removed == 1 || el_removed == 0);
      }
      else {
        /* forest_new consists only of complete trees. */
        T8_ASSERT (forest_new->incomplete_trees == 0);
        T8_ASSERT (forest_old->incomplete_trees == 0);
        /* If the levels differ, elem_new was refined or its family coarsened */
        if (level_old < level_new) {
          T8_ASSERT (level_new == level_old + 1);
          /* elem_old was refined */
          const t8_locidx_t family_size = scheme->element_get_num_children (tree_class, elem_old);
          const int refine = 1;
          replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, family_size, ielem_new);
          /* Advance to the next element */
          ielem_new += family_size;
          ielem_old++;
        }
        else if (level_old > level_new) {
          T8_ASSERT (level_new == level_old - 1);
          /* elem_old was coarsened */
          const t8_locidx_t family_size = scheme->element_get_num_children (tree_class, elem_new);
          const int refine = -1;
          replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, family_size, ielem_old, 1, ielem_new);
          /* Advance to the next element */
          ielem_new++;
          ielem_old += family_size;
        }
        else {
          /* elem_new = elem_old */
          T8_ASSERT (scheme->element_is_equal (tree_class, elem_new, elem_old));
          const int refine = 0;
          replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, 1, ielem_new);
          /* Advance to the next element */
          ielem_new++;
          ielem_old++;
        }
      }
    } /* element loop */
    T8_ASSERT (ielem_new == elems_per_tree_new);
    if (forest_new->incomplete_trees) {
      for (; ielem_old < elems_per_tree_old; ielem_old++) {
        /* remaining elements in old tree got removed */
        const int refine = -2;
        replace_fn (forest_old, forest_new, itree, tree_class, scheme, refine, 1, ielem_old, 0, -1);
      }
    }
    else {
      T8_ASSERT (ielem_old == elems_per_tree_old);
    }
  } /* tree loop */
  t8_global_productionf ("Done t8_forest_iterate_replace\n");
}

void
t8_forest_search_partition ([[maybe_unused]] const t8_forest_t forest,
                            [[maybe_unused]] t8_forest_partition_search_fn search_fn,
                            [[maybe_unused]] t8_forest_partition_query_fn query_fn,
                            [[maybe_unused]] sc_array_t *queries)
{
  SC_ABORT ("not implemented yet");
}

T8_EXTERN_C_END ();

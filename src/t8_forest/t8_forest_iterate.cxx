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

#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

typedef struct
{
  const t8_eclass_scheme_c *ts;
  int level;
  int num_children;
} t8_forest_child_type_query_t;

/* This is the function that we call in sc_split_array to determine for an
 * element E that is a descendant of an element e, of which of e's children, E is a descendant. */
static size_t
t8_forest_determine_child_type (sc_array_t *leaf_elements, size_t index, void *data)
{
  t8_forest_child_type_query_t *query_data = (t8_forest_child_type_query_t *) data;
  t8_element_t *element;

  /* Get a pointer to the element */
  element = (t8_element_t *) t8_sc_array_index_locidx (leaf_elements, index);
  T8_ASSERT (query_data->level < query_data->ts->t8_element_level (element));
  /* Compute the element's ancestor id at the stored level and return it as the element's type */
  return query_data->ts->t8_element_ancestor_id (element, query_data->level + 1);
}

void
t8_forest_split_array (const t8_element_t *element, t8_element_array_t *leaf_elements, size_t *offsets)
{
  sc_array_t offset_view;
  t8_forest_child_type_query_t query_data;

  const t8_eclass_scheme_c *ts = t8_element_array_get_scheme (leaf_elements);
  /* Store the number of children and the level of element */
  query_data.num_children = ts->t8_element_num_children (element);
  query_data.level = ts->t8_element_level (element);
  query_data.ts = ts;

  sc_array_t *element_array = t8_element_array_get_array_mutable (leaf_elements);
  /* Split the elements array according to the elements' ancestor id at
   * the given level. In other words for each child C of element, find
   * the indices i, j such that all descendants of C are
   * elements[i], ..., elements[j-1]
   */
  sc_array_init_data (&offset_view, offsets, sizeof (size_t), query_data.num_children + 1);
  sc_array_split (element_array, &offset_view, query_data.num_children, t8_forest_determine_child_type,
                  (void *) &query_data);
}

void
t8_forest_iterate_faces (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                         t8_element_array_t *leaf_elements, void *user_data, t8_locidx_t tree_lindex_of_first_leaf,
                         t8_forest_iterate_face_fn callback)
{
  t8_eclass_scheme_c *ts;
  t8_eclass_t eclass;
  t8_element_t **face_children;
  int child_face, num_face_children, iface;
  int *child_indices;
  size_t *split_offsets, indexa, indexb, elem_count;
  t8_element_array_t face_child_leaves;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));

  elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }
  eclass = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, eclass);

  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element
     * and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);
    T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
    if (ts->t8_element_equal (element, leaf)) {
      /* The element is the leaf, we are at the last stage of the recursion
       * and can call the callback. */
      (void) callback (forest, ltreeid, leaf, face, user_data, tree_lindex_of_first_leaf);
      return;
    }
  }
#ifdef T8_ENABLE_DEBUG
  /* Check whether element has greater level than the first leaf */
  const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);
  T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
  T8_ASSERT (ts->t8_element_level (element) < ts->t8_element_level (leaf));
#endif

  /* Call the callback function element, we pass -index - 1 as index to indicate
   * element is not a leaf, if it returns true, we continue with the top-down recursion */
  if (callback (forest, ltreeid, element, face, user_data, -tree_lindex_of_first_leaf - 1)) {
    /* Enter the recursion */
    /* We compute all face children of E, compute their leaf arrays and call iterate_faces */
    /* allocate the memory to store the face children */
    num_face_children = ts->t8_element_num_face_children (element, face);
    face_children = T8_ALLOC (t8_element_t *, num_face_children);
    ts->t8_element_new (num_face_children, face_children);
    /* Memory for the child indices of the face children */
    child_indices = T8_ALLOC (int, num_face_children);
    /* Memory for the indices that split the leaf_elements array */
    split_offsets = T8_ALLOC (size_t, ts->t8_element_num_children (element) + 1);
    /* Compute the face children */
    ts->t8_element_children_at_face (element, face, face_children, num_face_children, child_indices);
    /* Split the leaves array in portions belonging to the children of element */
    t8_forest_split_array (element, leaf_elements, split_offsets);
    for (iface = 0; iface < num_face_children; iface++) {
      /* Check if there are any leaf elements for this face child */
      indexa = split_offsets[child_indices[iface]];     /* first leaf of this face child */
      indexb = split_offsets[child_indices[iface] + 1]; /* first leaf of next child */
      if (indexa < indexb) {
        /* There exist leaves of this face child in leaf_elements,
         * we construct an array of these leaves */
        t8_element_array_init_view (&face_child_leaves, leaf_elements, indexa, indexb - indexa);
        /* Compute the corresponding face number of this face child */
        child_face = ts->t8_element_face_child_face (element, face, iface);
        /* Enter the recursion */
        t8_forest_iterate_faces (forest, ltreeid, face_children[iface], child_face, &face_child_leaves, user_data,
                                 indexa + tree_lindex_of_first_leaf, callback);
      }
    }
    /* clean-up */
    ts->t8_element_destroy (num_face_children, face_children);
    T8_FREE (face_children);
    T8_FREE (child_indices);
    T8_FREE (split_offsets);
  }
}

/* The recursion that is called from t8_forest_search_tree
 * Input is an element and an array of all leaf elements of this element.
 * The callback function is called on element and if it returns true,
 * the search continues with the children of the element.
 * Additionally a query function and a set of queries can be given.
 * In this case the recursion stops when either the search_fn function
 * returns false or the query_fn function returns false for all active queries.
 * (Thus, if there are no active queries left, the recursion also stops.)
 * A query is active for an element if the query_fn callback returned true
 * for the parent element.
 * If the callback function (search_fn) returns false for an element,
 * the query function is not called for this element.
 */
static void
t8_forest_search_recursion (t8_forest_t forest, const t8_locidx_t ltreeid, t8_element_t *element,
                            const t8_eclass_scheme_c *ts, t8_element_array_t *leaf_elements,
                            const t8_locidx_t tree_lindex_of_first_leaf, t8_forest_search_query_fn search_fn,
                            t8_forest_search_query_fn query_fn, sc_array_t *queries, sc_array_t *active_queries)
{
  /* Assertions to check for necessary requirements */
  /* The forest must be committed */
  T8_ASSERT (t8_forest_is_committed (forest));
  /* The tree must be local */
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));
  /* If we have queries, we also must have a query function */
  T8_ASSERT ((queries == NULL) == (query_fn == NULL));

  const size_t elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }
  const size_t num_active = queries == NULL ? 0 : active_queries->elem_count;
  if (queries != NULL && num_active == 0) {
    /* There are no queries left. We stop the recursion */
    return;
  }

  int is_leaf = 0;
  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);

    SC_CHECK_ABORT (ts->t8_element_level (element) <= ts->t8_element_level (leaf),
                    "Search: element level greater than leaf level\n");
    if (ts->t8_element_level (element) == ts->t8_element_level (leaf)) {
      T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
      T8_ASSERT (ts->t8_element_equal (element, leaf));
      /* The element is the leaf */
      is_leaf = 1;
    }
  }
  /* Call the callback function for the element */
  const int ret
    = search_fn (forest, ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf, NULL, NULL, NULL, 0);

  if (!ret) {
    /* The function returned false. We abort the recursion */
    return;
  }

  /* Check the queries.
   * If the current element is not a leaf, we store the queries that
   * return true in order to pass them on to the children of the element. */
  sc_array_t *new_active_queries = NULL;
  if (num_active > 0) {
    if (!is_leaf) {
      /* Initialize the new active query array */
      new_active_queries = sc_array_new (sizeof (size_t));
    }
    int *active_queries_matches = T8_ALLOC (int, num_active);
    T8_ASSERT (query_fn != NULL);
    query_fn (forest, ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf, queries, active_queries,
              active_queries_matches, num_active);

    for (size_t iactive = 0; iactive < num_active; iactive++) {
      if (!is_leaf && active_queries_matches[iactive]) {
        size_t query_index = *(size_t *) sc_array_index (active_queries, iactive);
        *(size_t *) sc_array_push (new_active_queries) = query_index;
      }
    }
    T8_FREE (active_queries_matches);
  }

  if (is_leaf) {
    /* The element was a leaf. We abort the recursion. */
    return;
  }

  if (num_active > 0 && new_active_queries->elem_count == 0) {
    /* No queries returned true for this element. We abort the recursion */
    sc_array_destroy (new_active_queries);
    return;
  }

  /* Enter the recursion (the element is definitely not a leaf at this point) */
  /* We compute all children of E, compute their leaf arrays and call search_recursion */
  /* allocate the memory to store the children */
  const int num_children = ts->t8_element_num_children (element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  ts->t8_element_new (num_children, children);
  /* Memory for the indices that split the leaf_elements array */
  size_t *split_offsets = T8_ALLOC (size_t, num_children + 1);
  /* Compute the children */
  ts->t8_element_children (element, num_children, children);
  /* Split the leaves array in portions belonging to the children of element */
  t8_forest_split_array (element, leaf_elements, split_offsets);
  for (int ichild = 0; ichild < num_children; ichild++) {
    /* Check if there are any leaf elements for this child */
    const size_t indexa = split_offsets[ichild];     /* first leaf of this child */
    const size_t indexb = split_offsets[ichild + 1]; /* first leaf of next child */
    if (indexa < indexb) {
      t8_element_array_t child_leaves;
      /* There exist leaves of this child in leaf_elements,
       * we construct an array of these leaves */
      t8_element_array_init_view (&child_leaves, leaf_elements, indexa, indexb - indexa);
      /* Enter the recursion */
      t8_forest_search_recursion (forest, ltreeid, children[ichild], ts, &child_leaves,
                                  indexa + tree_lindex_of_first_leaf, search_fn, query_fn, queries, new_active_queries);
    }
  }
  /* clean-up */
  ts->t8_element_destroy (num_children, children);
  T8_FREE (children);
  T8_FREE (split_offsets);
  if (num_active > 0) {
    sc_array_destroy (new_active_queries);
  }
}

/* Perform a top-down search in one tree of the forest */
static void
t8_forest_search_tree (t8_forest_t forest, t8_locidx_t ltreeid, t8_forest_search_query_fn search_fn,
                       t8_forest_search_query_fn query_fn, sc_array_t *queries, sc_array_t *active_queries)
{

  /* Get the element class, scheme and leaf elements of this tree */
  const t8_eclass_t eclass = t8_forest_get_eclass (forest, ltreeid);
  const t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest, eclass);
  t8_element_array_t *leaf_elements = t8_forest_tree_get_leaves (forest, ltreeid);

  /* assert for empty tree */
  T8_ASSERT (t8_element_array_get_count (leaf_elements) >= 0);
  /* Get the first and last leaf of this tree */
  const t8_element_t *first_el = t8_element_array_index_locidx (leaf_elements, 0);
  const t8_element_t *last_el
    = t8_element_array_index_locidx (leaf_elements, t8_element_array_get_count (leaf_elements) - 1);
  /* Compute their nearest common ancestor */
  t8_element_t *nca;
  ts->t8_element_new (1, &nca);
  ts->t8_element_nca (first_el, last_el, nca);

  /* Start the top-down search */
  t8_forest_search_recursion (forest, ltreeid, nca, ts, leaf_elements, 0, search_fn, query_fn, queries, active_queries);

  ts->t8_element_destroy (1, &nca);
}

void
t8_forest_search (t8_forest_t forest, t8_forest_search_query_fn search_fn, t8_forest_search_query_fn query_fn,
                  sc_array_t *queries)
{
  /* If we have queries build a list of all active queries,
   * thus all queries in the array */
  sc_array_t *active_queries = NULL;
  if (queries != NULL) {
    const size_t num_queries = queries->elem_count;
    /* build an array and write 0, 1, 2, 3,... into it */
    active_queries = sc_array_new_count (sizeof (size_t), num_queries);
    for (size_t iquery = 0; iquery < num_queries; ++iquery) {
      *(size_t *) sc_array_index (active_queries, iquery) = iquery;
    }
  }

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    t8_forest_search_tree (forest, itree, search_fn, query_fn, queries, active_queries);
  }

  if (active_queries != NULL) {
    sc_array_destroy (active_queries);
  }
}

void
t8_forest_iterate_replace (t8_forest_t forest_new, t8_forest_t forest_old, t8_forest_replace_t replace_fn)
{
  t8_global_productionf ("Into t8_forest_iterate_replace\n");
  T8_ASSERT (t8_forest_is_committed (forest_old));
  T8_ASSERT (t8_forest_is_committed (forest_new));

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest_new);
  T8_ASSERT (num_local_trees == t8_forest_get_num_local_trees (forest_old));

  for (t8_locidx_t itree = 0; itree < num_local_trees; itree++) {
    /* Loop over the trees */
    /* Get the number of elements of this tree in old and new forest */
    const t8_locidx_t elems_per_tree_new = t8_forest_get_tree_num_elements (forest_new, itree);
    const t8_locidx_t elems_per_tree_old = t8_forest_get_tree_num_elements (forest_old, itree);
    /* Get the eclass and scheme of the tree */
    t8_eclass_t eclass = t8_forest_get_tree_class (forest_new, itree);
    T8_ASSERT (eclass == t8_forest_get_tree_class (forest_old, itree));
    t8_eclass_scheme_c *ts = t8_forest_get_eclass_scheme (forest_new, eclass);
    T8_ASSERT (ts == t8_forest_get_eclass_scheme (forest_new, eclass));

    t8_locidx_t ielem_new = 0;
    t8_locidx_t ielem_old = 0;
    while (ielem_new < elems_per_tree_new) {
      /* Iterate over the elements */
      T8_ASSERT (ielem_new < elems_per_tree_new);
      T8_ASSERT (ielem_old < elems_per_tree_old);

      /* Get pointers to the elements */
      const t8_element_t *elem_new = t8_forest_get_element_in_tree (forest_new, itree, ielem_new);
      const t8_element_t *elem_old = t8_forest_get_element_in_tree (forest_old, itree, ielem_old);

      /* Get the levels of these elements */
      const int level_new = ts->t8_element_level (elem_new);
      const int level_old = ts->t8_element_level (elem_old);

      if (forest_new->incomplete_trees) {
        /* If el_removed is 1, the element in forest_new has been removed.
         * It is assumed that no element was removed. */
        int el_removed = 0;
        if (level_old < level_new) {
          /* elem_old got refined or removed */
          t8_element_t *elem_parent;
          ts->t8_element_new (1, &elem_parent);
          ts->t8_element_parent (elem_new, elem_parent);
          if (ts->t8_element_equal (elem_old, elem_parent)) {
            /* elem_old got refined */
            T8_ASSERT (level_new == level_old + 1);
            const t8_locidx_t family_size = ts->t8_element_num_children (elem_old);
#if T8_DEBUG
            /* Check if family of new refined elements is complete */
            T8_ASSERT (ielem_new + family_size <= elems_per_tree_new);
            for (t8_locidx_t ielem = 1; ielem < family_size; ielem++) {
              const t8_element_t *elem_new_debug = t8_forest_get_element_in_tree (forest_new, itree, ielem_new + ielem);
              ts->t8_element_parent (elem_new_debug, elem_parent);
              SC_CHECK_ABORT (ts->t8_element_equal (elem_old, elem_parent), "Family is not complete.");
            }
#endif
            ts->t8_element_destroy (1, &elem_parent);
            const int refine = 1;
            replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, family_size, ielem_new);
            /* Advance to the next element */
            ielem_new += family_size;
            ielem_old++;
          }
          else {
            /* elem_old got removed */
            el_removed = 1;
            ts->t8_element_destroy (1, &elem_parent);
          }
        }
        else if (level_old > level_new) {
          /* elem_old got coarsened or removed */
          t8_element_t *elem_parent;
          ts->t8_element_new (1, &elem_parent);
          ts->t8_element_parent (elem_old, elem_parent);
          if (ts->t8_element_equal (elem_new, elem_parent)) {
            /* elem_old got coarsened */
            T8_ASSERT (level_new == level_old - 1);
            /* Get size of family of old forest */
            int family_size = 1;
            for (t8_locidx_t ielem = 1;
                 ielem < ts->t8_element_num_children (elem_new) && ielem + ielem_old < elems_per_tree_old; ielem++) {
              elem_old = t8_forest_get_element_in_tree (forest_old, itree, ielem_old + ielem);
              ts->t8_element_parent (elem_old, elem_parent);
              if (ts->t8_element_equal (elem_new, elem_parent)) {
                family_size++;
              }
            }
            T8_ASSERT (family_size <= ts->t8_element_num_children (elem_new));
#if T8_DEBUG
            /* Check whether elem_old is the first element of the family */
            for (t8_locidx_t ielem = 1; ielem < ts->t8_element_num_children (elem_old) && ielem_old - ielem >= 0;
                 ielem++) {
              const t8_element_t *elem_old_debug = t8_forest_get_element_in_tree (forest_old, itree, ielem_old - ielem);
              ts->t8_element_parent (elem_old_debug, elem_parent);
              SC_CHECK_ABORT (!ts->t8_element_equal (elem_new, elem_parent),
                              "elem_old is not the first of the family.");
            }
#endif
            ts->t8_element_destroy (1, &elem_parent);
            const int refine = -1;
            replace_fn (forest_old, forest_new, itree, ts, refine, family_size, ielem_old, 1, ielem_new);
            /* Advance to the next element */
            ielem_new++;
            ielem_old += family_size;
          }
          else {
            /* elem_old got removed */
            el_removed = 1;
            ts->t8_element_destroy (1, &elem_parent);
          }
        }
        else {
          /* elem_old was untouched or got removed */
          if (ts->t8_element_equal (elem_new, elem_old)) {
            /* elem_new = elem_old */
            const int refine = 0;
            replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, 1, ielem_new);
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
          replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, 0, -1);
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
          const t8_locidx_t family_size = ts->t8_element_num_children (elem_old);
          const int refine = 1;
          replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, family_size, ielem_new);
          /* Advance to the next element */
          ielem_new += family_size;
          ielem_old++;
        }
        else if (level_old > level_new) {
          T8_ASSERT (level_new == level_old - 1);
          /* elem_old was coarsened */
          const t8_locidx_t family_size = ts->t8_element_num_children (elem_new);
          const int refine = -1;
          replace_fn (forest_old, forest_new, itree, ts, refine, family_size, ielem_old, 1, ielem_new);
          /* Advance to the next element */
          ielem_new++;
          ielem_old += family_size;
        }
        else {
          /* elem_new = elem_old */
          T8_ASSERT (ts->t8_element_equal (elem_new, elem_old));
          const int refine = 0;
          replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, 1, ielem_new);
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
        replace_fn (forest_old, forest_new, itree, ts, refine, 1, ielem_old, 0, -1);
      }
    }
    else {
      T8_ASSERT (ielem_old == elems_per_tree_old);
    }
  } /* tree loop */
  t8_global_productionf ("Done t8_forest_iterate_replace\n");
}

T8_EXTERN_C_END ();

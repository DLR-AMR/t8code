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

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/**
 * This struct stores query data about the elements within the forest.
 */
using t8_forest_child_type_query_t = struct
{
  const t8_scheme *scheme; /**< The scheme. */
  t8_eclass_t tree_class;  /**< The tree class. */
  int level;               /**< The refinement level. */
  int num_children;        /**< The number of children. */
};

/** 
 * This is the function that we call in sc_split_array to determine for an
 * element E that is a descendant of an element e, of which of e's children, E is a descendant. 
 * 
 * \param[in]   leaf_elements   The larray of eaf elements.
 * \param[in]   index           The local id of the element within the leaf_elements array.
 * \param[in]   data            The query data.
 * 
 * \return The element's ancestor id at the stored level is returned as the element's type.
 */
static size_t
t8_forest_determine_child_type (sc_array_t *leaf_elements, size_t index, void *data)
{
  t8_forest_child_type_query_t *query_data = (t8_forest_child_type_query_t *) data;
  t8_element_t *element;

  /* Get a pointer to the element */
  element = (t8_element_t *) t8_sc_array_index_locidx (leaf_elements, index);
  T8_ASSERT (query_data->level < query_data->scheme->element_get_level (query_data->tree_class, element));
  /* Compute the element's ancestor id at the stored level and return it as the element's type */
  return query_data->scheme->element_get_ancestor_id (query_data->tree_class, element, query_data->level + 1);
}

void
t8_forest_split_array (const t8_element_t *element, const t8_element_array_t *leaf_elements, size_t *offsets)
{
  sc_array_t offset_view;
  t8_forest_child_type_query_t query_data;

  const t8_scheme *scheme = t8_element_array_get_scheme (leaf_elements);
  const t8_eclass_t tree_class = t8_element_array_get_tree_class (leaf_elements);
  /* Store the number of children and the level of element */
  query_data.num_children = scheme->element_get_num_children (tree_class, element);
  query_data.level = scheme->element_get_level (tree_class, element);
  query_data.tree_class = tree_class;
  query_data.scheme = scheme;

  const sc_array_t *element_array = t8_element_array_get_array (leaf_elements);
  /* Split the elements array according to the elements' ancestor id at
   * the given level. In other words for each child C of element, find
   * the indices i, j such that all descendants of C are
   * elements[i], ..., elements[j-1]
   */
  sc_array_init_data (&offset_view, offsets, sizeof (size_t), query_data.num_children + 1);
  // Unfortunately, we have to cast away constness to pass to sc_array_split
  sc_array_split ((sc_array_t *) element_array, &offset_view, query_data.num_children, t8_forest_determine_child_type,
                  (void *) &query_data);
}

void
t8_forest_iterate_faces (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                         const int face, const t8_element_array_t *leaf_elements, void *user_data,
                         const t8_locidx_t tree_lindex_of_first_leaf, const t8_forest_iterate_face_fn callback)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
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

  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element
     * and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);
    T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
    if (scheme->element_is_equal (eclass, element, leaf)) {
      /* The element is the leaf, we are at the last stage of the recursion
       * and can call the callback. */
      (void) callback (forest, ltreeid, leaf, face, user_data, tree_lindex_of_first_leaf);
      return;
    }
  }
#if T8_ENABLE_DEBUG
  /* Check whether element has greater level than the first leaf */
  const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);
  T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
  T8_ASSERT (scheme->element_get_level (eclass, element) < scheme->element_get_level (eclass, leaf));
#endif

  /* Call the callback function element, we pass -index - 1 as index to indicate
   * element is not a leaf, if it returns true, we continue with the top-down recursion */
  if (callback (forest, ltreeid, element, face, user_data, -tree_lindex_of_first_leaf - 1)) {
    /* Enter the recursion */
    /* We compute all face children of E, compute their leaf arrays and call iterate_faces */
    /* allocate the memory to store the face children */
    num_face_children = scheme->element_get_num_face_children (eclass, element, face);
    face_children = T8_ALLOC (t8_element_t *, num_face_children);
    scheme->element_new (eclass, num_face_children, face_children);
    /* Memory for the child indices of the face children */
    child_indices = T8_ALLOC (int, num_face_children);
    /* Memory for the indices that split the leaf_elements array */
    split_offsets = T8_ALLOC (size_t, scheme->element_get_num_children (eclass, element) + 1);
    /* Compute the face children */
    scheme->element_get_children_at_face (eclass, element, face, face_children, num_face_children, child_indices);
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
        child_face = scheme->element_face_get_child_face (eclass, element, face, iface);
        /* Enter the recursion */
        t8_forest_iterate_faces (forest, ltreeid, face_children[iface], child_face, &face_child_leaves, user_data,
                                 indexa + tree_lindex_of_first_leaf, callback);
      }
    }
    /* clean-up */
    scheme->element_destroy (eclass, num_face_children, face_children);
    T8_FREE (face_children);
    T8_FREE (child_indices);
    T8_FREE (split_offsets);
  }
}

/**
 * The recursion that is called from t8_forest_search_tree.
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
 * 
 * \param[in] forest                    A committed forest.
 * \param[in] ltreeid                   The local index of the tree.
 * \param[in] element                   The considered element.
 * \param[in] tree_class                The element class of the tree.
 * \param[in] leaf_elements             The array of leaf elements that are descendants of \a element. Sorted by linear index.
 * \param[in] tree_lindex_of_first_leaf Tree-local index of the first leaf.
 * \param[in] search_fn                 The search callback function, \see t8_forest_search_fn.
 * \param[in] query_fn                  The query function \see t8_forest_query_fn.
 * \param[in] queries                   The array of queries.
 * \param[in] active_queries            The active queries, given as an array of indices within the \a queries array.
 */
static void
t8_forest_search_recursion (t8_forest_t forest, const t8_locidx_t ltreeid, t8_element_t *element,
                            const t8_eclass_t tree_class, t8_element_array_t *leaf_elements,
                            const t8_locidx_t tree_lindex_of_first_leaf, t8_forest_search_fn search_fn,
                            t8_forest_query_fn query_fn, sc_array_t *queries, sc_array_t *active_queries)
{
  /* Assertions to check for necessary requirements */
  /* The forest must be committed */
  T8_ASSERT (t8_forest_is_committed (forest));
  /* The tree must be local */
  T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));
  /* If we have queries, we also must have a query function */
  T8_ASSERT ((queries == NULL) == (query_fn == NULL));

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const size_t elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leaves left, so we have nothing to do */
    return;
  }
  const size_t num_active = queries == nullptr ? 0 : active_queries->elem_count;
  if (queries != nullptr && num_active == 0) {
    /* There are no queries left. We stop the recursion */
    return;
  }

  int is_leaf = 0;
  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element and if so call the callback function */
    const t8_element_t *leaf = t8_element_array_index_locidx (leaf_elements, 0);

    SC_CHECK_ABORT (scheme->element_get_level (tree_class, element) <= scheme->element_get_level (tree_class, leaf),
                    "Search: element level greater than leaf level\n");
    if (scheme->element_get_level (tree_class, element) == scheme->element_get_level (tree_class, leaf)) {
      T8_ASSERT (t8_forest_element_is_leaf (forest, leaf, ltreeid));
      T8_ASSERT (scheme->element_is_equal (tree_class, element, leaf));
      /* The element is the leaf */
      is_leaf = 1;
    }
  }
  /* Call the callback function for the element */
  const int ret = search_fn (forest, ltreeid, element, is_leaf, leaf_elements, tree_lindex_of_first_leaf);

  if (!ret) {
    /* The function returned false. We abort the recursion */
    return;
  }

  /* Check the queries.
   * If the current element is not a leaf, we store the queries that
   * return true in order to pass them on to the children of the element. */
  sc_array_t *new_active_queries = nullptr;
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
  const int num_children = scheme->element_get_num_children (tree_class, element);
  t8_element_t **children = T8_ALLOC (t8_element_t *, num_children);
  scheme->element_new (tree_class, num_children, children);
  /* Memory for the indices that split the leaf_elements array */
  size_t *split_offsets = T8_ALLOC (size_t, num_children + 1);
  /* Compute the children */
  scheme->element_get_children (tree_class, element, num_children, children);
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
      t8_forest_search_recursion (forest, ltreeid, children[ichild], tree_class, &child_leaves,
                                  indexa + tree_lindex_of_first_leaf, search_fn, query_fn, queries, new_active_queries);
    }
  }
  /* clean-up */
  scheme->element_destroy (tree_class, num_children, children);
  T8_FREE (children);
  T8_FREE (split_offsets);
  if (num_active > 0) {
    sc_array_destroy (new_active_queries);
  }
}

/** 
 *  Perform a top-down search in one tree of the forest 
 *  
 * \param[in] forest          The forest.
 * \param[in] ltreeid         The local index of the tree to search in.
 * \param[in] search_fn       The callback function describing the search criterion.
 * \param[in] query_fn        The callback function used for queries.
 * \param[in] queries         The array of queries.
 * \param[in] active_queries  The active queries within \a queries.
*/
static void
t8_forest_search_tree (t8_forest_t forest, t8_locidx_t ltreeid, t8_forest_search_fn search_fn,
                       t8_forest_query_fn query_fn, sc_array_t *queries, sc_array_t *active_queries)
{

  /* Get the element class, scheme and leaf elements of this tree */
  const t8_eclass_t eclass = t8_forest_get_eclass (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_element_array_t *leaf_elements = t8_forest_tree_get_leaf_elements (forest, ltreeid);

  /* Get the first and last leaf of this tree */
  const t8_element_t *first_el = t8_element_array_index_locidx (leaf_elements, 0);
  const t8_element_t *last_el
    = t8_element_array_index_locidx (leaf_elements, t8_element_array_get_count (leaf_elements) - 1);
  /* Compute their nearest common ancestor */
  t8_element_t *nca;
  scheme->element_new (eclass, 1, &nca);
  scheme->element_get_nca (eclass, first_el, last_el, nca);

  /* Start the top-down search */
  t8_forest_search_recursion (forest, ltreeid, nca, eclass, leaf_elements, 0, search_fn, query_fn, queries,
                              active_queries);

  scheme->element_destroy (eclass, 1, &nca);
}

void
t8_forest_search (t8_forest_t forest, t8_forest_search_fn search_fn, t8_forest_query_fn query_fn, sc_array_t *queries)
{
  /* If we have queries build a list of all active queries,
   * thus all queries in the array */
  sc_array_t *active_queries = nullptr;
  if (queries != nullptr) {
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

  if (active_queries != nullptr) {
    sc_array_destroy (active_queries);
  }
}

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
    t8_eclass_t tree_class = t8_forest_get_tree_class (forest_new, itree);
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

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
#include <t8_forest.h>
#include <t8_element_cxx.hxx>


/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

typedef struct
{
  t8_eclass_scheme_c *ts;
  int                 level;
  int                 num_children;
} t8_forest_child_type_query_t;

/* This is the function that we call in sc_split_array to determine for an
 * element E that is a descendant of an element e, of which of e's children,
 * E is a descendant. */
static              size_t
t8_forest_determine_child_type (sc_array_t * leaf_elements,
                                size_t index, void *data)
{
  t8_forest_child_type_query_t *query_data =
    (t8_forest_child_type_query_t *) data;
  t8_element_t       *element;

  /* Get a pointer to the element */
  element = (t8_element_t *) t8_sc_array_index_locidx (leaf_elements, index);
  T8_ASSERT (query_data->level < query_data->ts->t8_element_level (element));
  /* Compute the element's ancestor id at the stored level and return it
   * as the element's type */
  return query_data->ts->t8_element_ancestor_id (element,
                                                 query_data->level + 1);
}

void
t8_forest_split_array (const t8_element_t * element,
                       t8_element_array_t * leaf_elements, size_t * offsets)
{
  sc_array_t          offset_view;
  sc_array_t         *element_array;
  t8_forest_child_type_query_t query_data;
  t8_eclass_scheme_c *ts;

  ts = t8_element_array_get_scheme (leaf_elements);
  /* Store the number of children and the level of element */
  query_data.num_children = ts->t8_element_num_children (element);
  query_data.level = ts->t8_element_level (element);
  query_data.ts = ts;

  element_array = t8_element_array_get_array (leaf_elements);
  /* Split the elements array according to the elements' ancestor id at
   * the given level. In other words for each child C of element, find
   * the indices i, j such that all descendants of C are
   * elements[i], ..., elements[j-1]
   */
  sc_array_init_data (&offset_view, offsets, sizeof (size_t),
                      query_data.num_children + 1);
  sc_array_split (element_array, &offset_view, query_data.num_children,
                  t8_forest_determine_child_type, (void *) &query_data);
}

void
t8_forest_iterate_faces (t8_forest_t forest, t8_locidx_t ltreeid,
                         const t8_element_t * element, int face,
                         t8_element_array_t * leaf_elements, void *user_data,
                         t8_locidx_t tree_lindex_of_first_leaf,
                         t8_forest_iterate_face_fn callback)
{
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass;
  t8_element_t       *leaf, **face_children;
  int                 child_face, num_face_children, iface;
  int                *child_indices;
  size_t             *split_offsets, indexa, indexb, elem_count;
  t8_element_array_t  face_child_leafs;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid
             && ltreeid < t8_forest_get_num_local_trees (forest));

  elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leafs left, so we have nothing to do */
    return;
  }
  eclass = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, eclass);

  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element
     * and if so call the callback function */
    leaf = t8_element_array_index_locidx (leaf_elements, 0);
    if (!ts->t8_element_compare (element, leaf)) {
      /* The element is the leaf, we are at the last stage of the recursion
       * and can call the callback. */
      (void) callback (forest, ltreeid, leaf, face, user_data,
                       tree_lindex_of_first_leaf);
      return;
    }
  }
#ifdef T8_ENABLE_DEBUG
  /* Check whether element has greater level than the first leaf */
  leaf = t8_element_array_index_locidx (leaf_elements, 0);
  T8_ASSERT (ts->t8_element_level (element) < ts->t8_element_level (leaf));
#endif

  /* Call the callback function element, we pass -index - 1 as index to indicate
   * element is not a leaf, if it returns true, we continue with the
   * top-down recursion */
  if (callback (forest, ltreeid, element, face, user_data,
                -tree_lindex_of_first_leaf - 1)) {
    /* Enter the recursion */
    /* We compute all face children of E, compute their leaf arrays and
     * call iterate_faces */
    /* allocate the memory to store the face children */
    num_face_children = ts->t8_element_num_face_children (element, face);
    face_children = T8_ALLOC (t8_element_t *, num_face_children);
    ts->t8_element_new (num_face_children, face_children);
    /* Memory for the child indices of the face children */
    child_indices = T8_ALLOC (int, num_face_children);
    /* Memory for the indices that split the leaf_elements array */
    split_offsets =
      T8_ALLOC (size_t, ts->t8_element_num_children (element) + 1);
    /* Compute the face children */
    ts->t8_element_children_at_face (element, face, face_children,
                                     num_face_children, child_indices);
    /* Split the leafs array in portions belonging to the children of element */
    t8_forest_split_array (element, leaf_elements, split_offsets);
    for (iface = 0; iface < num_face_children; iface++) {
      /* Check if there are any leaf elements for this face child */
      indexa = split_offsets[child_indices[iface]];     /* first leaf of this face child */
      indexb = split_offsets[child_indices[iface] + 1]; /* first leaf of next child */
      if (indexa < indexb) {
        /* There exist leafs of this face child in leaf_elements,
         * we construct an array of these leafs */
        t8_element_array_init_view (&face_child_leafs, leaf_elements, indexa,
                                    indexb - indexa);
        /* Compute the corresponding face number of this face child */
        child_face = ts->t8_element_face_child_face (element, face, iface);
        /* Enter the recursion */
        t8_forest_iterate_faces (forest, ltreeid, face_children[iface],
                                 child_face, &face_child_leafs, user_data,
                                 indexa + tree_lindex_of_first_leaf,
                                 callback);
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
 * the search continues with the children of the element. */
static void
t8_forest_search_recursion (t8_forest_t forest, t8_locidx_t ltreeid,
                            t8_eclass_t eclass, t8_element_t * element,
                            t8_eclass_scheme_c * ts,
                            t8_element_array_t * leaf_elements,
                            t8_locidx_t tree_lindex_of_first_leaf,
                            t8_forest_search_query_fn search_fn,
                            void *user_data)
{
  t8_element_t       *leaf, **children;
  int                 num_children, ichild;
  size_t             *split_offsets, indexa, indexb;
  t8_element_array_t  child_leafs;
  size_t              elem_count;
  int                 ret;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid
             && ltreeid < t8_forest_get_num_local_trees (forest));

  elem_count = t8_element_array_get_count (leaf_elements);
  if (elem_count == 0) {
    /* There are no leafs left, so we have nothing to do */
    return;
  }

  if (elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element
     * and if so call the callback function */
    leaf = t8_element_array_index_locidx (leaf_elements, 0);

    SC_CHECK_ABORT (ts->t8_element_level (element) <=
                    ts->t8_element_level (leaf),
                    "Search: element level greater than leaf level\n");
    if (ts->t8_element_level (element) == ts->t8_element_level (leaf)) {
      T8_ASSERT (!ts->t8_element_compare (element, leaf));
      /* The element is the leaf, we are at the last stage of the recursion
       * and can call the callback. */
      (void) search_fn (forest, ltreeid, leaf, leaf_elements, user_data,
                        tree_lindex_of_first_leaf);
      return;
    }
  }
  /* Call the callback function for the element, we pass -index -1 as index to indicate
   * element is not a leaf */
  ret = search_fn (forest, ltreeid, element, leaf_elements, user_data,
                   -tree_lindex_of_first_leaf - 1);

  if (ret) {
    /* Enter the recursion */
    /* We compute all children of E, compute their leaf arrays and
     * call search_recursion */
    /* allocate the memory to store the children */
    num_children = ts->t8_element_num_children (element);
    children = T8_ALLOC (t8_element_t *, num_children);
    ts->t8_element_new (num_children, children);
    /* Memory for the indices that split the leaf_elements array */
    split_offsets = T8_ALLOC (size_t, num_children + 1);
    /* Compute the children */
    ts->t8_element_children (element, num_children, children);
    /* Split the leafs array in portions belonging to the children of element */
    t8_forest_split_array (element, leaf_elements, split_offsets);
    for (ichild = 0; ichild < num_children; ichild++) {
      /* Check if there are any leaf elements for this child */
      indexa = split_offsets[ichild];   /* first leaf of this child */
      indexb = split_offsets[ichild + 1];       /* first leaf of next child */
      if (indexa < indexb) {
        /* There exist leafs of this child in leaf_elements,
         * we construct an array of these leafs */
        t8_element_array_init_view (&child_leafs, leaf_elements, indexa,
                                    indexb - indexa);
        /* Enter the recursion */
        t8_forest_search_recursion (forest, ltreeid, eclass, children[ichild],
                                    ts, &child_leafs,
                                    indexa + tree_lindex_of_first_leaf,
                                    search_fn, user_data);
      }
    }
    /* clean-up */
    ts->t8_element_destroy (num_children, children);
    T8_FREE (children);
    T8_FREE (split_offsets);
  }
}

/* Perform a top-down search in one tree of the forest */
static void
t8_forest_search_tree (t8_forest_t forest, t8_locidx_t ltreeid,
                       t8_forest_search_query_fn search_fn, void *user_data)
{
  t8_eclass_t         eclass;
  t8_eclass_scheme_c *ts;
  t8_element_t       *nca, *first_el, *last_el;
  t8_element_array_t *leaf_elements;

  /* Get the element class, scheme and leaf elements of this tree */
  eclass = t8_forest_get_eclass (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, eclass);
  leaf_elements = t8_forest_tree_get_leafs (forest, ltreeid);

  /* assert for empty tree */
  T8_ASSERT (t8_element_array_get_count (leaf_elements) >= 0);
  /* Get the first and last leaf of this tree */
  first_el = t8_element_array_index_locidx (leaf_elements, 0);
  last_el =
    t8_element_array_index_locidx (leaf_elements,
                                   t8_element_array_get_count (leaf_elements)
                                   - 1);
  /* Compute their nearest common ancestor */
  ts->t8_element_new (1, &nca);
  ts->t8_element_nca (first_el, last_el, nca);
  /* Start the top-down search */
  t8_forest_search_recursion (forest, ltreeid, eclass, nca, ts, leaf_elements,
                              0, search_fn, user_data);
}

void
t8_forest_search (t8_forest_t forest, t8_forest_search_query_fn search_fn,
                  void *user_data)
{
  t8_locidx_t         num_local_trees, itree;

  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_local_trees; itree++) {
    t8_forest_search_tree (forest, itree, search_fn, user_data);
  }
}

void
t8_forest_iterate_replace (t8_forest_t forest_new,
                           t8_forest_t forest_old,
                           t8_forest_replace_t replace_fn)
{
  t8_locidx_t         ielem_new, ielem_old, elems_per_tree_old,
    elems_per_tree_new;
  t8_locidx_t         itree, num_local_trees;
  t8_locidx_t         family_size;
  t8_element_t       *elem_new, *elem_old;
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass;
  int                 level_new, level_old;

  t8_global_productionf ("Into t8_forest_iterate_replace\n");
  T8_ASSERT (t8_forest_is_committed (forest_old));
  T8_ASSERT (t8_forest_is_committed (forest_new));

  num_local_trees = t8_forest_get_num_local_trees (forest_new);
  T8_ASSERT (num_local_trees == t8_forest_get_num_local_trees (forest_old));

  for (itree = 0; itree < num_local_trees; itree++) {
    /* Loop over the trees */
    /* Get the number of elements of this tree in old and new forest */
    elems_per_tree_new = t8_forest_get_tree_num_elements (forest_new, itree);
    elems_per_tree_old = t8_forest_get_tree_num_elements (forest_old, itree);
    /* Get the eclass and scheme of the tree */
    eclass = t8_forest_get_tree_class (forest_new, itree);
    T8_ASSERT (eclass == t8_forest_get_tree_class (forest_old, itree));
    ts = t8_forest_get_eclass_scheme (forest_new, eclass);
    T8_ASSERT (ts == t8_forest_get_eclass_scheme (forest_new, eclass));
    for (ielem_new = 0, ielem_old = 0; ielem_new < elems_per_tree_new
         || ielem_old < elems_per_tree_old;) {
      /* Iterate over the elements */
      /* Get pointers to the elements */
      elem_new = t8_forest_get_element_in_tree (forest_new, itree, ielem_new);
      elem_old = t8_forest_get_element_in_tree (forest_old, itree, ielem_old);
      /* Get the levels of these elements */
      level_new = ts->t8_element_level (elem_new);
      level_old = ts->t8_element_level (elem_old);
      /* If the levels differ, elem_new was refined or its family coarsened */
      if (level_old < level_new) {
        T8_ASSERT (level_new == level_old + 1);
        /* elem_old was refined */
        family_size = ts->t8_element_num_children (elem_old);
        replace_fn (forest_old, forest_new, itree, ts, 1, ielem_old,
                    family_size, ielem_new);
        /* Advance to the next element */
        ielem_new += family_size;
        ielem_old++;
      }
      else if (level_old > level_new) {
        T8_ASSERT (level_new == level_old - 1);
        /* elem_old was coarsened */
        family_size = ts->t8_element_num_children (elem_new);
        replace_fn (forest_old, forest_new, itree, ts, family_size, ielem_old,
                    1, ielem_new);
        /* Advance to the next element */
        ielem_new++;
        ielem_old += family_size;
      }
      else {
        /* elem_new = elem_old */
        T8_ASSERT (!ts->t8_element_compare (elem_new, elem_old));
        replace_fn (forest_old, forest_new, itree, ts, 1, ielem_old, 1,
                    ielem_new);
        /* Advance to the next element */
        ielem_new++;
        ielem_old++;
      }
    }                           /* element loop */
    T8_ASSERT (ielem_new ==
               t8_forest_get_tree_num_elements (forest_new, itree));
    T8_ASSERT (ielem_old ==
               t8_forest_get_tree_num_elements (forest_old, itree));
  }                             /* tree loop */
  t8_global_productionf ("Done t8_forest_iterate_replace\n");
}

T8_EXTERN_C_END ();

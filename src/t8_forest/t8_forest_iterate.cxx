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

#include <t8_default/t8_dtri.h>

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
t8_forest_determine_child_type (sc_array_t * leaf_elements, size_t index,
                                void *data)
{
  t8_forest_child_type_query_t *query_data =
    (t8_forest_child_type_query_t *) data;
  t8_element_t       *element;

  /* Get a pointer to the element */
  element = (t8_element_t *) sc_array_index (leaf_elements, index);
  T8_ASSERT (query_data->level < query_data->ts->t8_element_level (element));
  /* Compute the element's ancestor id at the stored level and return it
   * as the element's type */
  return query_data->ts->t8_element_ancestor_id (element,
                                                 query_data->level + 1);
}

void
t8_forest_split_array (const t8_element_t * element,
                       sc_array_t * leaf_elements, t8_eclass_scheme_c * ts,
                       size_t * offsets)
{
  sc_array_t          offset_view;
  t8_forest_child_type_query_t query_data;

  /* Store the number of children and the level of element */
  query_data.num_children = ts->t8_element_num_children (element);
  query_data.level = ts->t8_element_level (element);
  query_data.ts = ts;

  /* Split the elements array according to the elements' ancestor id at
   * the given level. In other words for each child C of element, find
   * the indices i, j such that all descendants of C are
   * elements[i], ..., elements[j-1]
   */
  sc_array_init_data (&offset_view, offsets, sizeof (size_t),
                      query_data.num_children + 1);
  sc_array_split (leaf_elements, &offset_view, query_data.num_children,
                  t8_forest_determine_child_type, (void *) &query_data);
}

void
t8_forest_iterate_faces (t8_forest_t forest, t8_locidx_t ltreeid,
                         const t8_element_t * element, int face,
                         sc_array_t * leaf_elements, void *user_data,
                         t8_locidx_t tree_lindex_of_first_leaf,
                         t8_forest_iterate_face_fn callback)
{
  t8_eclass_scheme_c *ts;
  t8_eclass_t         eclass;
  t8_element_t       *leaf, **face_children;
  int                 child_face, num_face_children, iface;
  int                *child_indices;
  size_t             *split_offsets, indexa, indexb;
  sc_array_t          face_child_leafs;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (0 <= ltreeid
             && ltreeid < t8_forest_get_num_local_trees (forest));

  if (leaf_elements->elem_count == 0) {
    /* There are no leafs left, so we have nothing to do */
    return;
  }
  eclass = t8_forest_get_tree_class (forest, ltreeid);
  ts = t8_forest_get_eclass_scheme (forest, eclass);

  if (leaf_elements->elem_count == 1) {
    /* There is only one leaf left, we check whether it is the same as element
     * and if so call the callback function */
    leaf = (t8_element_t *) sc_array_index (leaf_elements, 0);
    if (!ts->t8_element_compare (element, leaf)) {
      /* The element is the leaf, we are at the last stage of the recursion
       * and can call the callback. */
      callback (forest, ltreeid, leaf, face, user_data);
      return;
    }
  }
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
  t8_debugf ("[H] Split array for element with level %i\n",
             ts->t8_element_level (element));
  t8_forest_split_array (element, leaf_elements, ts, split_offsets);
  for (iface = 0; iface < num_face_children; iface++) {
    /* Check if there are any leaf elements for this face child */
    indexa = split_offsets[child_indices[iface]];       /* first leaf of this face child */
    indexb = split_offsets[child_indices[iface] + 1];   /* first leaf of next child */
    t8_debugf ("[H] On face child %i, child_index %i, ina %i inb %i\n",
               iface, child_indices[iface], indexa, indexb);
    if (indexa < indexb) {
      /* There exist leafs of this face child in leaf_elements,
       * we construct an array of these leafs */
      sc_array_init_view (&face_child_leafs, leaf_elements, indexa,
                          indexb - indexa);
      /* Compute the corresponding face number of this face child */
      child_face = ts->t8_element_face_child_face (element, face, iface);
      t8_debugf ("[H] Call rec with face child %i level %i\n",
                 iface, ts->t8_element_level (face_children[iface]));
      /* Enter the recursion */
      t8_forest_iterate_faces (forest, ltreeid, face_children[iface],
                               child_face, &face_child_leafs, user_data,
                               callback);
    }
  }
  /* clean-up */
  ts->t8_element_destroy (num_face_children, face_children);
  T8_FREE (face_children);
  T8_FREE (child_indices);
  T8_FREE (split_offsets);
}

T8_EXTERN_C_END ();

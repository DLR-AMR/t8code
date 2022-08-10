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

#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest.h>
#include <t8_data/t8_containers.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

#if T8_ENABLE_DEBUG
/** Return zero if the first \a num_elemets in \a elements are not a (sub)family.
 * \param [in] tscheme       The element scheme for current local tree 
 *                           where the elements are from.
 * \param [in] elements      The elements array.
 * \param [in] num_elements  The first \a num_elements to be checked in \a elements.
 * \note If the first element has level 0, the return is 0.
 * \note This test does not compare with the elements before and after the elements
 *       given by \a elements in the current forest. A non-zero return is therefore
 *       not valid. 
 */
static int
t8_forest_is_family_callback (t8_eclass_scheme_c *ts,
                              const int num_elements, 
                              t8_element_t **elements) {
  t8_element_t       *element_parent;
  t8_element_t       *element_parent_compare;
  int                 iter;
  int                 num_siblings;
  int                 level;

  for (iter = 0; iter < num_elements; iter++) {
    T8_ASSERT (ts->t8_element_is_valid (elements[iter]));
  }

  level = ts->t8_element_level(elements[0]);
  if (level == 0) {
    return 0;
  }

  ts->t8_element_new (1, &element_parent_compare);
  ts->t8_element_new (1, &element_parent);
  ts->t8_element_parent (elements[0], element_parent);
  num_siblings = ts->t8_element_num_siblings (elements[0]);

  for (iter = 0; iter < num_elements; iter++) {
    ts->t8_element_parent (elements[iter], element_parent_compare);
    if (ts->t8_element_compare (element_parent, element_parent_compare)) {
      ts->t8_element_destroy (1, &element_parent);
      ts->t8_element_destroy (1, &element_parent_compare);
      return 0;
    }
  }

  if (num_elements < num_siblings) {
    for (iter = num_elements; iter < num_elements; iter++) {
      if (elements[iter] != NULL) {
        T8_ASSERT (ts->t8_element_is_valid (elements[iter]));
        ts->t8_element_parent (elements[iter], element_parent_compare);
        if (!ts->t8_element_compare (element_parent, element_parent_compare)) {
          ts->t8_element_destroy (1, &element_parent);
          ts->t8_element_destroy (1, &element_parent_compare);
          return 0;
        }
      }
    }
  }

  ts->t8_element_destroy (1, &element_parent);
  ts->t8_element_destroy (1, &element_parent_compare);
  return 1;
}
#endif

/** Return the index of the first family member of a given family in an array of elements.
 * \param [in] forest        The forest
 * \param [in] ts            The element scheme for elements in \a telements.
 * \param [in] telements     The array of newly created (adapted) elements.
 * \param [in] telements_pos The index of an element in \a telement
 *                           array which could be coarsened recursively.
 * \return                   The index of the first family member whose family is
 *                           defined by \a telements_pos in \a telement.                    
 * \note The element with index \a telements_pos must be the last child in its family.
 *       \see t8_forest_adapt_coarsen_recursive.
 * \note If the element with index \a telements_pos in \a telement can not be coarsened
 *       recursively, return INT32_MIN.
 */
static t8_locidx_t
t8_forest_pos (t8_forest_t forest,
               t8_eclass_scheme_c *ts,
               t8_element_array_t *telements,
               t8_locidx_t telements_pos)
{
  t8_element_t       *element; /* element with index telements_pos */
  t8_element_t       *element_compare;
  t8_element_t       *element_parent;
  t8_element_t       *element_parent_compare;
  int                 level_current;
  int                 level_compare;
  t8_locidx_t         child_id;
  t8_locidx_t         pos;
  t8_locidx_t         el_iter; /* Loop running variable */
  t8_locidx_t         elements_in_array;
  t8_locidx_t         num_siblings;

  pos = telements_pos;
  elements_in_array = t8_element_array_get_count (telements);
  T8_ASSERT (0 <= pos && pos < elements_in_array);
  element = t8_element_array_index_locidx (telements, pos);
  element_compare = t8_element_array_index_locidx (telements, pos);
  level_current = ts->t8_element_level (element);
  num_siblings = ts->t8_element_num_siblings (element);
  child_id = ts->t8_element_child_id (element);
  
  /* Elements whose family consist of exactly one element do not 
   * get coarsened recursively.
   * If child_id is not last, elements cannot be coarsened recursively.
   * Elements with level 0 cannot be further coarsened. */
  if (!(child_id > 0 && child_id == num_siblings - 1) || level_current == 0) {
    return INT32_MIN;
  }
  /* If the forest is complete, the family is also complete. 
   * Thus, the index of the first member can be determined. */
  if (!forest->is_incomplete) {
    return pos - num_siblings - 1;
  }
  
  T8_ASSERT (child_id > 0 && child_id == num_siblings - 1);
  T8_ASSERT (level_current > 0);

  ts->t8_element_new (1, &element_parent_compare);
  ts->t8_element_new (1, &element_parent);
  /* Get parent of a family member by coarsening last member. */
  ts->t8_element_parent (element, element_parent);

  /* Loop backward over all possible familie members until we hit an 
   * element that is not part of the family or we have reached the 
   * maximum number of member. */
  for (el_iter = 1; el_iter < num_siblings 
                    && el_iter < elements_in_array; el_iter++) {
    pos = telements_pos - el_iter;
    T8_ASSERT (0 <= pos && pos < elements_in_array);
    element_compare = t8_element_array_index_locidx (telements, pos);
    /* Quick check by level. Not mandatory. */
    level_compare = ts->t8_element_level (element_compare);
    if (level_current != level_compare) {
      break;
    }
    ts->t8_element_parent (element_compare, element_parent_compare);
    if (ts->t8_element_compare (element_parent, element_parent_compare)) {
      break;
    }
  }

  /* If the current family is smaller than possible, check that the first 
   * element along the space-filling-curve next to the family is overlapped 
   * when family is coarsened. */
  if (el_iter < num_siblings && el_iter < elements_in_array) {
    level_compare = ts->t8_element_level (element_compare);
    /* Only elements with higher level then level of elements in family, can get 
     * potentially be overlapped. */
    if (level_compare > level_current) {
      /* Compare ancestors */
      ts->t8_element_nca (element, element_compare, element_parent_compare);
      level_compare = ts->t8_element_level (element_parent_compare);
      T8_ASSERT (level_compare <= level_current-1);
      if (level_compare == level_current-1) {
        pos = INT32_MIN; /* No recursion coarsening */
      }
    }
    pos++;
  }     

  /* clean up */
  ts->t8_element_destroy (1, &element_parent);
  ts->t8_element_destroy (1, &element_parent_compare);

#if T8_ENABLE_MPI
  /* The first element on process rank must have child_id 0, otherwise other 
   * family members could be on process rank-1. */
  if (pos == 0 && forest->mpirank > 0) {
    element = t8_element_array_index_locidx (telements, pos);
    child_id = ts->t8_element_child_id (element);
    if (child_id > 0) {
      return INT32_MIN;
    }
  }
#endif

  return pos;
}

/** Check the lastly inserted elements of an array for recursive coarsening.
 * The last inserted element must be the last element of a family.
 * \param [in] forest  The new forest currently in construction.
 * \param [in] ltreeid The current local tree.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] ts      The scheme for this local tree.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in] el_coarsen the index of the first element in \a telement
 *                        array which could be coarsened recursively.
 * \param [in,out] el_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Buffer space to store a family of elements.
 */
static void
t8_forest_adapt_coarsen_recursive (t8_forest_t forest,
                                   t8_locidx_t ltreeid,
                                   t8_locidx_t lelement_id,
                                   t8_eclass_scheme_c *ts,
                                   t8_element_array_t *telements,
                                   t8_locidx_t el_coarsen,
                                   t8_locidx_t *el_inserted,
                                   t8_element_t **el_buffer)
{
  t8_element_t       *element;
  t8_element_t      **fam;
  t8_locidx_t         pos;
  t8_locidx_t         i; /* Loop running variable */
  t8_locidx_t         elements_in_array;
  t8_locidx_t         num_siblings;
  t8_locidx_t         num_elements_to_adapt_callback;
  int                 is_family;
  /* el_inserted is the index of the last element in telements plus one.
   * el_coarsen is the index of the first element which could possibly
   * be coarsened. */
  element = t8_element_array_index_locidx (telements, *el_inserted - 1);
  elements_in_array = t8_element_array_get_count (telements);
  num_siblings = ts->t8_element_num_siblings (element);
  fam = el_buffer;
  is_family = 1;
  T8_ASSERT (*el_inserted == (t8_locidx_t) elements_in_array);
  T8_ASSERT (ts->t8_element_level (element) > 0);
  T8_ASSERT (el_coarsen >= 0);
  pos = t8_forest_pos (forest, ts, telements, *el_inserted - 1);

  while (is_family && pos >= el_coarsen && pos < elements_in_array) {
    /* Get all elements at indices pos, pos + 1, ... ,pos + num_siblings - 1 */
#if T8_ENABLE_DEBUG
    for (i = 0; i < num_siblings; i++) {
      fam[i] = NULL;
    }
#endif
    for (i = 0; i < num_siblings && pos + i < elements_in_array; i++) {
      fam[i] = t8_element_array_index_locidx (telements, pos + i);
    }

    if (forest->is_incomplete) {
      /* We will pass a (in)complete family to the adapt callback */
      num_elements_to_adapt_callback = *el_inserted - pos;
      T8_ASSERT (0 < num_elements_to_adapt_callback);
      T8_ASSERT (num_elements_to_adapt_callback <= num_siblings);
    }
    else if (i == num_siblings && ts->t8_element_is_family (fam)) {
      /* We will pass a full family to the adapt callback */
      num_elements_to_adapt_callback = num_siblings;
    }
    else {
      is_family = 0;
      num_elements_to_adapt_callback = 1;
    }
    T8_ASSERT (forest->is_incomplete ? (!is_family || t8_forest_is_family_callback 
                                         (ts, num_elements_to_adapt_callback, fam))
              : (!is_family || ts->t8_element_is_family (fam)));
    if (is_family
        && forest->set_adapt_fn (forest, forest->set_from, ltreeid,
                                 lelement_id, ts, is_family, 
                                 num_elements_to_adapt_callback, fam) == -1) {
      /* Coarsen the element */
      *el_inserted -= num_elements_to_adapt_callback - 1;
      /* remove num_elements_to_adapt_callback - 1 elements from the array */
      T8_ASSERT ((size_t) elements_in_array == t8_element_array_get_count (telements));
      T8_ASSERT (ts->t8_element_level 
                      (t8_element_array_index_locidx (telements, pos)) > 0);
      ts->t8_element_parent (fam[0], fam[0]);
      /*Shorten the array by the number of siblings of the fine element */
      elements_in_array -= num_elements_to_adapt_callback - 1;
      num_siblings = ts->t8_element_num_siblings (fam[0]);
      t8_element_array_resize (telements, elements_in_array);
      /* Set element to the new constructed parent. Since resizing the array
       * may change the position in memory, we have to do it after resizing. */
      T8_ASSERT (*el_inserted-1 == pos);
      pos = t8_forest_pos (forest, ts, telements, pos);
    }
    else {
      /* If the elements are no family or
       * the family is not to be coarsened we abort the coarsening process */
      is_family = 0;
    }
  } /* End while loop */
}

/** Check the lastly inserted element of an array for recursive refining.
 * \param [in] forest  The new forest currently in construction.
 * \param [in] ltreeid The current local tree.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] ts      The scheme for this local tree.
 * \param [in] elem_list Helper list to temporarily insert the newly refined elements.
 *                       These will eventually get copied to \a telements.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in,out] num_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Enough buffer space to store all children of the lastly created element.
 */
static void
t8_forest_adapt_refine_recursive (t8_forest_t forest, t8_locidx_t ltreeid,
                                  t8_locidx_t lelement_id,
                                  t8_eclass_scheme_c *ts,
                                  sc_list_t * elem_list,
                                  t8_element_array_t *telements,
                                  t8_locidx_t *num_inserted,
                                  t8_element_t **el_buffer)
{
  t8_element_t       *insert_el;
  int                 num_children;
  int                 ci;
  int                 refine;

  if (elem_list->elem_count <= 0) {
    return;
  }

  while (elem_list->elem_count > 0) {
    /* Until the list is empty we
     * - remove the first element from the list.
     * - Check whether it should get refined or removed
     * - If refined, we add all its children to the list
     * - If removed, we just remove it from the list
     * - Otherwise, we add the element to the array of new elements
     */
    el_buffer[0] = (t8_element_t *) sc_list_pop (elem_list);
    num_children = ts->t8_element_num_children (el_buffer[0]);
    refine = forest->set_adapt_fn (forest, forest->set_from, ltreeid, lelement_id,
                                   ts, 0, 1, el_buffer);
    T8_ASSERT (refine != -1);
    if (refine == 1) {
      /* The element should be refined */
      if (ts->t8_element_level (el_buffer[0]) < forest->maxlevel) {
        /* only refine, if we do not exceed the maximum allowed level */
        /* Create the children and add them to the list */
        ts->t8_element_new (num_children - 1, el_buffer + 1);
        ts->t8_element_children (el_buffer[0], num_children, el_buffer);
        for (ci = num_children - 1; ci >= 0; ci--) {
          (void) sc_list_prepend (elem_list, el_buffer[ci]);
        }
      }
    }
    else if (refine == -2) {
      /* This element should get removed,
       * we just remove it from the buffer */
      ts->t8_element_destroy (1, el_buffer);
    }
    else {
      T8_ASSERT (refine == 0);
      /* This element should not get refined,
       * we remove it from the buffer and add it to the array of new elements. */
      insert_el = t8_element_array_push (telements);
      ts->t8_element_copy (el_buffer[0], insert_el);
      ts->t8_element_destroy (1, el_buffer);
      (*num_inserted)++;
    }
  } /* End while loop */
}

/* TODO: optimize this when we own forest_from */
void
t8_forest_adapt (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  t8_eclass_scheme_c *tscheme;
  t8_element_array_t *telements;
  t8_element_array_t *telements_from;
  t8_element_t      **elements;
  t8_element_t      **elements_from;
  t8_locidx_t         ltree_id;
  t8_locidx_t         num_trees;
  t8_locidx_t         num_el_from;
  t8_locidx_t         el_considered;
  t8_locidx_t         el_inserted;
  t8_locidx_t         el_coarsen;
  t8_locidx_t         el_offset;
  t8_tree_t           tree;
  t8_tree_t           tree_from;
  sc_list_t          *refine_list = NULL;       /* This is only needed when we adapt recursively */
  size_t              num_children; 
  size_t              num_siblings; 
  size_t              zz; 
  size_t              curr_size_elements_from;
  size_t              curr_size_elements;
  size_t              num_elements_to_adapt_callback;
  int                 ci;
  int                 refine;
  int                 is_family;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->set_adapt_recursive != -1);

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start adadpt %f %f\n", sc_MPI_Wtime (),
                           forest->profile->adapt_runtime);
  }

  forest_from = forest->set_from;
  t8_global_productionf ("Into t8_forest_adapt from %lld total elements\n",
                         (long long) forest_from->global_num_elements);

  /* TODO: Allocate memory for the trees of forest.
   * Will we do this here or in an extra function? */
  T8_ASSERT (forest->trees->elem_count == forest_from->trees->elem_count);

  if (forest->set_adapt_recursive) {
    refine_list = sc_list_new (NULL);
  }
  forest->local_num_elements = 0;
  el_offset = 0;
  num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over the trees and build the new element arrays for each one. */
  for (ltree_id = 0; ltree_id < num_trees; ltree_id++) {
    /* Get the new and old tree and the new and old element arrays */
    tree = t8_forest_get_tree (forest, ltree_id);
    tree_from = t8_forest_get_tree (forest_from, ltree_id);
    telements = &tree->elements;
    telements_from = &tree_from->elements;
    const t8_element_t *first_element_from = t8_element_array_index_locidx
      (telements_from, 0);
    /* Number of elements in the old tree */
    num_el_from = (t8_locidx_t) t8_element_array_get_count (telements_from);
    T8_ASSERT (num_el_from ==
               t8_forest_get_tree_num_elements (forest_from, ltree_id));
    /* Get the element scheme for this tree */
    tscheme = t8_forest_get_eclass_scheme (forest_from, tree->eclass);
    /* Index of the element we currently consider for refinement/coarsening. */
    el_considered = 0;
    /* Index into the newly inserted elements */
    el_inserted = 0;
    /* el_coarsen is the index of the first element in the new element
     * array which could be coarsened recursively. */
    el_coarsen = 0;
    num_children = tscheme->t8_element_num_children (first_element_from);
    curr_size_elements = num_children;
    curr_size_elements_from =
      tscheme->t8_element_num_siblings (first_element_from);
    /* Buffer for a family of new elements */
    elements           = T8_ALLOC (t8_element_t *, num_children);
    /* Buffer for a family of old elements */
    elements_from      = T8_ALLOC (t8_element_t *, curr_size_elements_from);
    /* We now iterate over all elements in this tree and check them for refinement/coarsening. */
    while (el_considered < num_el_from) {
      /* Load the current element and at most num_siblings-1 many others into
       * the elements_from buffer. Stop when we are certain that they cannot from
       * a family.
       * At the end is_family will be true, if these elements form a family.
       */

      num_siblings =
        tscheme->t8_element_num_siblings (t8_element_array_index_locidx
                                          (telements_from, el_considered));

      if (num_siblings > curr_size_elements_from) {
        /* Enlarge the elements_from buffer if required */
        elements_from =
          T8_REALLOC (elements_from, t8_element_t *, num_siblings);
        curr_size_elements_from = num_siblings;
      }
#if T8_ENABLE_DEBUG
      for (zz = 0; zz < (unsigned int) num_siblings; zz++) {
        elements_from[zz] = NULL;
      }
#endif
      for (zz = 0; zz < (unsigned int) num_siblings &&
           el_considered + (t8_locidx_t) zz < num_el_from; zz++) {
        elements_from[zz] = t8_element_array_index_locidx (telements_from,
                                                           el_considered +
                                                           zz);
        /* This is a quick check whether we build up a family here and could
         * abort early if not.
         * If the child id of the current element is not zz, then it cannot
         * be part of a family (Since we can only have a family if child ids
         * are 0, 1, 2, ... zz, ... num_siblings-1).
         * This check is however not sufficient - therefore, we call is_family later. */
        if (!forest_from->is_incomplete &&
            (size_t) tscheme->t8_element_child_id (elements_from[zz]) != zz) {
          break;
        }
      }

      /* We assume that the elements do not form a family.
       * So we will only pass the first element to the adapt callback. */
      is_family = 0;
      num_elements_to_adapt_callback = 1;
      if (forest_from->is_incomplete) {
        is_family = 
          t8_forest_is_incomplete_family (forest_from, ltree_id, el_considered,
                                          tscheme, elements_from, zz);
        if (is_family > 0) {
          /* We will pass a (in)complete family to the adapt callback */
          num_elements_to_adapt_callback = is_family;
          is_family = 1;
        }
      }
      else if (zz == num_siblings
               && tscheme->t8_element_is_family (elements_from)) {
        /* We will pass a full family to the adapt callback */
        is_family = 1;
        num_elements_to_adapt_callback = num_siblings;
      }
      T8_ASSERT (num_elements_to_adapt_callback <= num_siblings);
      T8_ASSERT (forest_from->is_incomplete ? (!is_family || t8_forest_is_family_callback 
                                (tscheme, num_elements_to_adapt_callback, elements_from))
                : (!is_family || tscheme->t8_element_is_family (elements_from)));
      /* Pass the element, or the family to the adapt callback.
       * The output will be  1 if the element should be refined
       *                     0 if the element should remain as is
       *                    -1 if we passed a family and it should get coarsened
       *                    -2 if the element should be removed.
       */
      refine =
        forest->set_adapt_fn (forest, forest->set_from, ltree_id,
                              el_considered, tscheme, is_family,
                              num_elements_to_adapt_callback, elements_from);

      T8_ASSERT (is_family || refine != -1);
      if (refine > 0 && tscheme->t8_element_level (elements_from[0]) >=
          forest->maxlevel) {
        /* Only refine an element if it does not exceed the maximum level */
        refine = 0;
      }
      if (refine == 1) {
        /* The first element is to be refined */
        num_children = tscheme->t8_element_num_children (elements_from[0]);
        if (num_children > curr_size_elements) {
          elements = T8_REALLOC (elements, t8_element_t *, num_children);
          curr_size_elements = num_children;
        }
        if (forest->set_adapt_recursive) {
          /* Create the children of this element */
          tscheme->t8_element_new (num_children, elements);
          tscheme->t8_element_children (elements_from[0], num_children,
                                        elements);
          for (ci = num_children - 1; ci >= 0; ci--) {
            /* Prepend the children to the refine_list.
             * These should now be the only elements in the list.
             */
            (void) sc_list_prepend (refine_list, elements[ci]);
          }
          /* We now recursively check the newly created elements for refinement. */
          t8_forest_adapt_refine_recursive (forest, ltree_id, el_considered,
                                            tscheme, refine_list, telements,
                                            &el_inserted, elements);
          /* el_coarsen is the index of the first element in the new element
           * array which could be coarsened recursively.
           * We can set this here to the next element after the current family, 
           * since a family that emerges from a refinement will never be coarsened */
          el_coarsen = el_inserted + num_children;
        }
        else {
          (void) t8_element_array_push_count (telements, num_children);
          for (zz = 0; zz < num_children; zz++) {
            elements[zz] =
              t8_element_array_index_locidx (telements, el_inserted + zz);
          }
          tscheme->t8_element_children (elements_from[0], num_children,
                                        elements);
          el_inserted += num_children;
        }
        el_considered++;
      }
      else if (refine == -1) {
        /* The elements form a family and are to be coarsened. */
        /* Make room for one more new element. */
        elements[0] = t8_element_array_push (telements);
        /* Compute the parent of the current family.
         * This parent is now inserted in telements. */
        T8_ASSERT (tscheme->t8_element_level (elements_from[0]) > 0);
        tscheme->t8_element_parent (elements_from[0], elements[0]);
        /*num_siblings is now equivalent to the number of children of elements[0],
         * as num_siblings is always associated with elements_from*/
        num_children = num_siblings;
        el_inserted++;
        if (num_children > curr_size_elements) {
          elements = T8_REALLOC (elements, t8_element_t *, num_children);
          curr_size_elements = num_children;
        }
        if (forest->set_adapt_recursive) {
          /* Adaptation is recursive.
           * We check whether the just generated parent is the last in its
           * family (and not the only one).
           * If so, we check this family for recursive coarsening. */
          const int           child_id =
            tscheme->t8_element_child_id (elements[0]);
          if (child_id > 0 && (size_t) child_id == num_children - 1) {
            t8_forest_adapt_coarsen_recursive (forest, ltree_id,
                                               el_considered, tscheme,
                                               telements, el_coarsen,
                                               &el_inserted, elements);
          }
        }
        el_considered += num_elements_to_adapt_callback;    
      }
      else if (refine == 0) {
        /* The considered elements are neither to be coarsened nor is the first
         * one to be refined.
         * We copy the element to the new element array. */
        elements[0] = t8_element_array_push (telements);
        tscheme->t8_element_copy (elements_from[0], elements[0]);
        el_inserted++;
        if (forest->set_adapt_recursive) {
          /* Adaptation is recursive.
           * If adaptation is recursive and this was the last element in its family
           * (and not the only one), we need to check for recursive coarsening. */
          const int           child_id =
            tscheme->t8_element_child_id (elements[0]);
          if (child_id > 0 && (size_t) child_id == num_children - 1) {
            t8_forest_adapt_coarsen_recursive (forest, ltree_id, el_considered,
                                               tscheme, telements, el_coarsen,
                                               &el_inserted, elements);
          }
        } 
        el_considered++;
      }
      else {
        T8_ASSERT (refine == -2);
        /* Remove the element */
        el_considered++;
      }
    } /* End element loop */

    /* Check that if we had recursive adaptation, the refine list is now empty. */
    T8_ASSERT (!forest->set_adapt_recursive || refine_list->elem_count == 0);
#if 1
    /* Empty trees, even subtrees on processes, lead to problems. 
     * When all elements have been removed from new tree, insert the last 
     * element from old tree (tree_from). */
    if (!el_inserted) {
      T8_ASSERT (refine == -2 || (refine == 1 && forest->set_adapt_recursive));
      T8_ASSERT (!(t8_locidx_t) t8_element_array_get_count (telements));
      /* We copy the last element to the new element array. */
      elements[0] = t8_element_array_push (telements);
      tscheme->t8_element_copy (elements_from[0], elements[0]);
      el_inserted++;
    }
    /* Set the new element offset of this tree */
    tree->elements_offset = el_offset;
    el_offset += el_inserted;
    /* Add to the new number of local elements. */
    forest->local_num_elements += el_inserted;
    /* Possibly shrink the telements array to the correct size */
    t8_element_array_resize (telements, el_inserted);
#else
    /* Empty trees lead to problems. When all elements have been removed from 
     * new tree, remove this tree. */
    if (!el_inserted) {
      T8_ASSERT (refine == -2 || (refine == 1 && forest->set_adapt_recursive));
      T8_ASSERT(!(t8_locidx_t) t8_element_array_get_count (telements));
      /* remove empty tree from forests */
      t8_forest_remove_tree (forest, ltree_id);
      /* TODO: We need extra variable num_trees_removed otherwise the tree loop 
       * is lost in ballance. */
      num_trees_removed++;
      num_trees--;
      ltree_id--;
    }
    else {
      /* Set the new element offset of this tree */
      tree->elements_offset = el_offset;
      el_offset += el_inserted;
      /* Add to the new number of local elements. */
      forest->local_num_elements += el_inserted;
      /* Possibly shrink the telements array to the correct size */
      t8_element_array_resize (telements, el_inserted);
    }
#endif
    /* clean up */
    T8_FREE (elements);
    T8_FREE (elements_from);
  } /* End tree loop */
  if (forest->set_adapt_recursive) {
    /* clean up */
    sc_list_destroy (refine_list);
  }

  /* We now adapted all local trees */
  /* Compute the new global number of elements */
  t8_forest_comm_global_num_elements (forest);
  t8_global_productionf ("Done t8_forest_adapt with %lld total elements\n",
                         (long long) forest->global_num_elements);
  
  SC_CHECK_ABORT (forest->global_num_elements > 0, "Forest is empty.");

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime += sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End adadpt %f %f\n", sc_MPI_Wtime (),
                           forest->profile->adapt_runtime);
  }

}

T8_EXTERN_C_END ();
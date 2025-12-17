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
/** \file t8_forest_adapt.cxx
 * Implements functions declared in \ref t8_forest_adapt.h.
 */

#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_data/t8_containers.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

#if T8_ENABLE_DEBUG
/** Return zero if the first \a num_elements in \a elements are not a (sub)family.
 * \param [in] scheme       The element scheme for current local tree 
 *                           where the elements are from.
 * \param [in] tree_class    The eclass of tree the elements are part of.
 * \param [in] elements      The elements array.
 * \param [in] num_elements  The first \a num_elements to be checked in \a elements.
 * \return                   0 if the first \a num_elements in \a elements are 
 *                           not a (sub)family.
 * \note If the first element has level 0, the return is 0.
 * \note This test does not compare with the elements before and after the elements
 *       given by \a elements in the current forest. A non-zero return is therefore
 *       not valid. 
 */
static int
t8_forest_is_family_callback (const t8_scheme *scheme, t8_eclass_t tree_class, const int num_elements,
                              t8_element_t **elements)
{

  for (int iter = 0; iter < num_elements; iter++) {
    T8_ASSERT (scheme->element_is_valid (tree_class, elements[iter]));
  }

  if (scheme->element_get_level (tree_class, elements[0]) == 0) {
    return 0;
  }

  t8_element_t *element_parent;
  t8_element_t *element_parent_compare;
  scheme->element_new (tree_class, 1, &element_parent_compare);
  scheme->element_new (tree_class, 1, &element_parent);
  scheme->element_get_parent (tree_class, elements[0], element_parent);

  for (int iter = 0; iter < num_elements; iter++) {
    scheme->element_get_parent (tree_class, elements[iter], element_parent_compare);
    if (!scheme->element_is_equal (tree_class, element_parent, element_parent_compare)) {
      scheme->element_destroy (tree_class, 1, &element_parent);
      scheme->element_destroy (tree_class, 1, &element_parent_compare);
      return 0;
    }
  }

  if (num_elements < scheme->element_get_num_siblings (tree_class, elements[0])) {
    for (int iter = num_elements; iter < num_elements; iter++) {
      scheme->element_get_parent (tree_class, elements[iter], element_parent_compare);
      if (scheme->element_is_equal (tree_class, element_parent, element_parent_compare)) {
        scheme->element_destroy (tree_class, 1, &element_parent);
        scheme->element_destroy (tree_class, 1, &element_parent_compare);
        return 0;
      }
    }
  }

  scheme->element_destroy (tree_class, 1, &element_parent);
  scheme->element_destroy (tree_class, 1, &element_parent_compare);
  return 1;
}
#endif

/** Return the index of the first family member of a given family in an array of elements.
 * \param [in] forest        The forest
 * \param [in] tree_class    The eclass of tree the elements are part of.
 * \param [in] scheme        The element scheme for elements in \a telements.
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
t8_forest_pos (t8_forest_t forest, t8_eclass_t tree_class, const t8_scheme *scheme, t8_element_array_t *telements,
               const t8_locidx_t telements_pos)
{
#if T8_ENABLE_DEBUG
  const t8_eclass_t tree_eclass_check = telements->tree_class;
  T8_ASSERT (tree_eclass_check == tree_class);
#endif

  const t8_locidx_t elements_in_array = t8_element_array_get_count (telements);
  T8_ASSERT (0 <= telements_pos && telements_pos < elements_in_array);
  const t8_element_t *element = t8_element_array_index_locidx (telements, telements_pos);
  const int level_current = scheme->element_get_level (tree_class, element);
  const int num_siblings = scheme->element_get_num_siblings (tree_class, element);

  {
    const int child_id = scheme->element_get_child_id (tree_class, element);
    /* Left if condition:
     * If child_id is not last, elements cannot be coarsened recursively.
     * But elements (vertex) whose family consist of exactly one element do 
     * also not get coarsened recursively.
     * Right if condition:
     * Elements with level 0 cannot be further coarsened. */
    if (!(child_id > 0 && child_id == num_siblings - 1) || level_current == 0) {
      return INT32_MIN;
    }
    T8_ASSERT (child_id > 0 && child_id == num_siblings - 1);
    T8_ASSERT (level_current > 0);
  }

  /* If the forest is complete, the family is also complete. 
   * Thus, the index of the first member can be determined. */
  if (!forest->incomplete_trees) {
    return telements_pos - (t8_locidx_t) num_siblings - 1;
  }

  t8_element_t *element_parent;
  t8_element_t *element_parent_compare;
  scheme->element_new (tree_class, 1, &element_parent_compare);
  scheme->element_new (tree_class, 1, &element_parent);
  /* Get parent of a family member by coarsening last member. */
  scheme->element_get_parent (tree_class, element, element_parent);

  /* Loop backward over all possible family members until we hit an 
   * element that is not part of the family or we have reached the 
   * maximum number of member. */
  t8_locidx_t el_iter; /* Loop running variable */
  t8_locidx_t pos = -1;
  const t8_element_t *element_compare;
  for (el_iter = 1; el_iter < (t8_locidx_t) num_siblings && el_iter < elements_in_array; el_iter++) {
    pos = telements_pos - el_iter;
    T8_ASSERT (0 <= pos && pos < elements_in_array);
    element_compare = t8_element_array_index_locidx (telements, pos);
    const int level_compare = scheme->element_get_level (tree_class, element_compare);
    /* By comparing the levels in advance we may be able to avoid
     * the more complex test with the parent element.*/
    if (level_current != level_compare) {
      break;
    }
    scheme->element_get_parent (tree_class, element_compare, element_parent_compare);
    if (!scheme->element_is_equal (tree_class, element_parent, element_parent_compare)) {
      break;
    }
  }

  /* If the current set of considered elements is smaller in size than a possible 
   * family, check if the first element along the space-filling-curve next to the
   * considered elements is overlapped when set is coarsened. */
  if (el_iter < (t8_locidx_t) num_siblings && el_iter < elements_in_array) {
    int level_compare = scheme->element_get_level (tree_class, element_compare);
    /* Only elements with higher level then level of elements in family, can get 
     * potentially be overlapped. */
    if (level_compare > level_current) {
      /* Compare ancestors */
      scheme->element_get_nca (tree_class, element, element_compare, element_parent_compare);
      level_compare = scheme->element_get_level (tree_class, element_parent_compare);
      T8_ASSERT (level_compare <= level_current - 1);
      if (level_compare == level_current - 1) {
        pos = INT32_MIN; /* No recursion coarsening */
      }
    }
    pos++;
  }
  /* clean up */
  scheme->element_destroy (tree_class, 1, &element_parent);
  scheme->element_destroy (tree_class, 1, &element_parent_compare);

#if T8_ENABLE_MPI
  /* The first element on process rank must have child_id 0, otherwise other 
   * family members could be on process rank-1. */
  if (pos == 0 && forest->mpirank > 0) {
    const t8_element_t *element_boarder = t8_element_array_index_locidx (telements, pos);
    const int child_id = scheme->element_get_child_id (tree_class, element_boarder);
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
 * \param [in] tree_class The eclass of tree \a ltreeid.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] scheme      The scheme for this local tree.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in] el_coarsen the index of the first element in \a telement
 *                        array which could be coarsened recursively.
 * \param [in,out] el_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Buffer space to store a family of elements.
 */
static void
t8_forest_adapt_coarsen_recursive (t8_forest_t forest, t8_locidx_t ltreeid, t8_eclass_t tree_class,
                                   t8_locidx_t lelement_id, const t8_scheme *scheme, t8_element_array_t *telements,
                                   t8_locidx_t el_coarsen, t8_locidx_t *el_inserted, t8_element_t **el_buffer)
{
  T8_ASSERT (el_coarsen >= 0);
  /* el_inserted is the index of the last element in telements plus one.
   * el_coarsen is the index of the first element which could possibly
   * be coarsened. */
  int num_siblings;
  {
    const t8_element_t *element = t8_element_array_index_locidx (telements, *el_inserted - 1);
    T8_ASSERT (scheme->element_get_level (tree_class, element) > 0);
    num_siblings = scheme->element_get_num_siblings (tree_class, element);
  }

  t8_locidx_t elements_in_array = t8_element_array_get_count (telements);
  T8_ASSERT (*el_inserted == (t8_locidx_t) elements_in_array);
  t8_element_t **fam = el_buffer;
  int is_family = 1;
  t8_locidx_t pos = t8_forest_pos (forest, tree_class, scheme, telements, *el_inserted - 1);

  while (is_family && pos >= el_coarsen && pos < elements_in_array) {
    t8_locidx_t ielement; /* Loop running variable */
                          /* Get all elements at indices pos, pos + 1, ... ,pos + num_siblings - 1 */
#if T8_ENABLE_DEBUG
    for (ielement = 0; ielement < (t8_locidx_t) num_siblings; ielement++) {
      fam[ielement] = NULL;
    }
#endif
    for (ielement = 0; ielement < (t8_locidx_t) num_siblings && pos + ielement < elements_in_array; ielement++) {
      /* TODO: In a future version, fam[ielement] should be const and we should call t8_element_array_index_locidx (the const version). */
      fam[ielement] = t8_element_array_index_locidx_mutable (telements, pos + ielement);
    }

    int num_elements_to_adapt_callback;
    if (forest->set_from->incomplete_trees) {
      /* We will pass a (in)complete family to the adapt callback */
      num_elements_to_adapt_callback = (int) (*el_inserted - pos);
      T8_ASSERT (0 < num_elements_to_adapt_callback);
      T8_ASSERT (num_elements_to_adapt_callback <= num_siblings);
    }
    else if (ielement == (t8_locidx_t) num_siblings && scheme->elements_are_family (tree_class, fam)) {
      /* We will pass a full family to the adapt callback */
      num_elements_to_adapt_callback = num_siblings;
    }
    else {
      is_family = 0;
      num_elements_to_adapt_callback = 1;
    }
#if T8_ENABLE_DEBUG
    /* If is_family is true, the set fam must be a family. */
    if (forest->set_from->incomplete_trees) {
      T8_ASSERT (!is_family || t8_forest_is_family_callback (scheme, tree_class, num_elements_to_adapt_callback, fam));
    }
    else {
      T8_ASSERT (forest->set_from->incomplete_trees == 0);
      T8_ASSERT (!is_family || scheme->elements_are_family (tree_class, fam));
    }
#endif
    if (is_family
        && forest->set_adapt_fn (forest->set_from, ltreeid, tree_class, lelement_id, scheme, is_family,
                                 num_elements_to_adapt_callback, fam, forest->user_data, forest->t8code_data)
             == -1) {
      /* Coarsen the element */
      *el_inserted -= (t8_locidx_t) (num_elements_to_adapt_callback - 1);
      /* remove num_elements_to_adapt_callback - 1 elements from the array */
      T8_ASSERT ((size_t) elements_in_array == t8_element_array_get_count (telements));
      T8_ASSERT (scheme->element_get_level (tree_class, t8_element_array_index_locidx (telements, pos)) > 0);
      scheme->element_get_parent (tree_class, fam[0], fam[0]);
      /*Shorten the array by the number of siblings of the fine element */
      elements_in_array -= (t8_locidx_t) num_elements_to_adapt_callback - 1;
      num_siblings = scheme->element_get_num_siblings (tree_class, fam[0]);
      t8_element_array_resize (telements, elements_in_array);
      /* Set element to the new constructed parent. Since resizing the array
       * may change the position in memory, we have to do it after resizing. */
      T8_ASSERT (*el_inserted - 1 == pos);
      pos = t8_forest_pos (forest, tree_class, scheme, telements, pos);
    }
    else {
      /* If the elements are no family or
       * the family is not to be coarsened we abort the coarsening process */
      is_family = 0;
    }
  } /* End while loop */
}

/** Check the lastly inserted element of an array for recursive refining or removing.
 * \param [in] forest  The new forest currently in construction.
 * \param [in] ltreeid The current local tree.
 * \param [in] tree_class The eclass of tree \a ltreeid.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] scheme      The scheme for this local tree.
 * \param [in] elem_list Helper list to temporarily insert the newly refined elements.
 *                       These will eventually get copied to \a telements.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in,out] num_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Enough buffer space to store all children of the lastly created element.
 * \param [in] element_removed Flag set to 1 if element was removed.
 */
static void
t8_forest_adapt_refine_recursive (t8_forest_t forest, t8_locidx_t ltreeid, t8_eclass_t tree_class,
                                  t8_locidx_t lelement_id, const t8_scheme *scheme, sc_list_t *elem_list,
                                  t8_element_array_t *telements, t8_locidx_t *num_inserted, t8_element_t **el_buffer,
                                  int *element_removed)
{
  while (elem_list->elem_count > 0) {
    /* Until the list is empty we
     * - remove the first element from the list.
     * - Check whether it should get refined or removed
     * - If refined, we add all its children to the list
     * - If removed, we just remove it from the list
     * - Otherwise, we add the element to the array of new elements
     */
    el_buffer[0] = (t8_element_t *) sc_list_pop (elem_list);
    const int num_children = scheme->element_get_num_children (tree_class, el_buffer[0]);
    const int is_family = 0;
    const int num_elements_to_adapt_callback = 1;
    const int refine
      = forest->set_adapt_fn (forest->set_from, ltreeid, tree_class, lelement_id, scheme, is_family,
                              num_elements_to_adapt_callback, el_buffer, forest->user_data, forest->t8code_data);
    T8_ASSERT (refine != -1);
    if (refine == 1) {
      /* The element should be refined */
      if (scheme->element_get_level (tree_class, el_buffer[0]) < forest->maxlevel
          && scheme->element_is_refinable (tree_class, el_buffer[0])) {
        /* only refine if element is refinable and if we do not exceed the maximum allowed level */
        /* Create the children and add them to the list */
        scheme->element_new (tree_class, num_children - 1, el_buffer + 1);
        scheme->element_get_children (tree_class, el_buffer[0], num_children, el_buffer);
        for (int ci = num_children - 1; ci >= 0; ci--) {
          (void) sc_list_prepend (elem_list, el_buffer[ci]);
        }
      }
    }
    else if (refine == -2) {
      /* This element should get removed,
       * we just remove it from the buffer */
      *element_removed = 1;
      scheme->element_destroy (tree_class, 1, el_buffer);
    }
    else {
      T8_ASSERT (refine == 0);
      /* This element should not get refined,
       * we remove it from the buffer and add it to the array of new elements. */
      t8_element_t *insert_el = t8_element_array_push (telements);
      scheme->element_copy (tree_class, el_buffer[0], insert_el);
      scheme->element_destroy (tree_class, 1, el_buffer);
      (*num_inserted)++;
    }
  } /* End while loop */
}

/* TODO: optimize this when we own forest_from */
void
t8_forest_adapt (t8_forest_t forest)
{
  t8_forest_t forest_from;
  t8_element_array_t *telements;
  t8_element_array_t *telements_from;
  t8_element_t **elements;
  t8_element_t **elements_from;
  t8_locidx_t ltree_id;
  t8_locidx_t num_trees;
  t8_locidx_t num_el_from;
  t8_locidx_t el_considered;
  t8_locidx_t el_inserted;
  t8_locidx_t el_coarsen;
  t8_locidx_t el_offset;
  t8_tree_t tree;
  t8_tree_t tree_from;
  sc_list_t *refine_list = NULL; /* This is only needed when we adapt recursively */
  int num_children;
  int num_siblings;
  int curr_size_elements_from;
  int num_elements_to_adapt_callback;
  int zz;
  int ci;
  int refine;
  int is_family;
  int element_removed = 0;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->set_adapt_recursive != -1);

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start adapt %f %f\n", sc_MPI_Wtime (), forest->profile->adapt_runtime);
  }

  forest_from = forest->set_from;
  t8_global_productionf ("Into t8_forest_adapt from %lld total elements\n",
                         (long long) forest_from->global_num_leaf_elements);

  T8_ASSERT (forest_from->incomplete_trees != -1);
  T8_ASSERT (forest->incomplete_trees == -1);
  /* TODO: Allocate memory for the trees of forest.
   * Will we do this here or in an extra function? */
  T8_ASSERT (forest->trees->elem_count == forest_from->trees->elem_count);

  if (forest->set_adapt_recursive) {
    refine_list = sc_list_new (NULL);
  }
  forest->local_num_leaf_elements = 0;
  el_offset = 0;
  num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over the trees and build the new element arrays for each one. */
  for (ltree_id = 0; ltree_id < num_trees; ltree_id++) {
    /* Get the new and old tree and the new and old element arrays */
    tree = t8_forest_get_tree (forest, ltree_id);
    tree_from = t8_forest_get_tree (forest_from, ltree_id);
    telements = &tree->leaf_elements;
    telements_from = &tree_from->leaf_elements;
    /* Number of elements in the old tree */
    num_el_from = (t8_locidx_t) t8_element_array_get_count (telements_from);
    T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
    /* Continue only if tree_from is not empty.
     * Otherwise there is nothing to adapt, since elements can't be inserted. */
    if (num_el_from > 0) {
      const t8_element_t *first_element_from = t8_element_array_index_locidx (telements_from, 0);
      /* Get the element scheme for this tree */
      const t8_scheme *scheme = t8_forest_get_scheme (forest_from);
      /* Index of the element we currently consider for refinement/coarsening. */
      el_considered = 0;
      /* Index into the newly inserted elements */
      el_inserted = 0;
      /* el_coarsen is the index of the first element in the new element
       * array which could be coarsened recursively. */
      el_coarsen = 0;
      num_children = scheme->get_max_num_children (tree->eclass);
      curr_size_elements_from = scheme->element_get_num_siblings (tree->eclass, first_element_from);
      /* Buffer for a family of new elements */
      elements = T8_ALLOC (t8_element_t *, num_children);
      /* Buffer for a family of old elements */
      elements_from = T8_ALLOC (t8_element_t *, curr_size_elements_from);
      /* We now iterate over all elements in this tree and check them for refinement/coarsening. */
      while (el_considered < num_el_from) {
        /* Load the current element and at most num_siblings-1 many others into
         * the elements_from buffer. Stop when we are certain that they cannot from
         * a family.
         * At the end is_family will be true, if these elements form a family.
         */

        num_siblings = scheme->element_get_num_siblings (tree->eclass,
                                                         t8_element_array_index_locidx (telements_from, el_considered));

        if (num_siblings > curr_size_elements_from) {
          /* Enlarge the elements_from buffer if required */
          elements_from = T8_REALLOC (elements_from, t8_element_t *, num_siblings);
          curr_size_elements_from = num_siblings;
        }
#if T8_ENABLE_DEBUG
        for (zz = 0; zz < num_siblings; zz++) {
          elements_from[zz] = NULL;
        }
#endif
        for (zz = 0; zz < num_siblings && el_considered + (t8_locidx_t) zz < num_el_from; zz++) {
          /* TODO: In a future version elements_from[zz] should be const and we should call t8_element_array_index_locidx (the const version). */
          elements_from[zz] = t8_element_array_index_locidx_mutable (telements_from, el_considered + (t8_locidx_t) zz);
          /* This is a quick check whether we build up a family here and could
           * abort early if not.
           * If the child id of the current element is not zz, then it cannot
           * be part of a family (Since we can only have a family if child ids
           * are 0, 1, 2, ... zz, ... num_siblings-1).
           * This check is however not sufficient - therefore, we call is_family later. */
          if (!forest_from->incomplete_trees && scheme->element_get_child_id (tree->eclass, elements_from[zz]) != zz) {
            break;
          }
        }

        /* We assume that the elements do not form a family.
         * So we will only pass the first element to the adapt callback. */
        is_family = 0;
        num_elements_to_adapt_callback = 1;
        if (forest_from->incomplete_trees) {
          is_family = t8_forest_is_incomplete_family (forest_from, ltree_id, el_considered, elements_from, zz);
          if (is_family > 0) {
            /* We will pass a (in)complete family to the adapt callback */
            num_elements_to_adapt_callback = is_family;
            is_family = 1;
          }
        }
        else if (zz == num_siblings && scheme->elements_are_family (tree->eclass, elements_from)) {
          /* We will pass a full family to the adapt callback */
          is_family = 1;
          num_elements_to_adapt_callback = num_siblings;
        }
        T8_ASSERT (num_elements_to_adapt_callback <= num_siblings);
#if T8_ENABLE_DEBUG
        if (forest_from->incomplete_trees) {
          T8_ASSERT (forest_from->incomplete_trees == 1);
          T8_ASSERT (
            !is_family
            || t8_forest_is_family_callback (scheme, tree->eclass, num_elements_to_adapt_callback, elements_from));
        }
        else {
          T8_ASSERT (forest_from->incomplete_trees == 0);
          T8_ASSERT (!is_family || scheme->elements_are_family (tree->eclass, elements_from));
        }
#endif
        /* Pass the element, or the family to the adapt callback.
         * The output will be  1 if the element should be refined
         *                     0 if the element should remain as is
         *                    -1 if we passed a family and it should get coarsened
         *                    -2 if the element should be removed.
         */
        refine = forest->set_adapt_fn (forest->set_from, ltree_id, tree->eclass, el_considered, scheme, is_family,
                                       num_elements_to_adapt_callback, elements_from, forest->user_data,
                                       forest->t8code_data);

        T8_ASSERT (is_family || refine != -1);
        if (refine > 0
            && (scheme->element_get_level (tree->eclass, elements_from[0]) >= forest->maxlevel
                || !scheme->element_is_refinable (tree->eclass, elements_from[0]))) {
          /* Only refine an element if it does not exceed the maximum level and if it is refinable */
          refine = 0;
        }
        if (refine == 1) {
          /* The first element is to be refined */
          num_children = scheme->element_get_num_children (tree->eclass, elements_from[0]);
          if (forest->set_adapt_recursive) {
            /* Create the children of this element */
            scheme->element_new (tree->eclass, num_children, elements);
            scheme->element_get_children (tree->eclass, elements_from[0], num_children, elements);
            for (ci = num_children - 1; ci >= 0; ci--) {
              /* Prepend the children to the refine_list.
               * These should now be the only elements in the list.
               */
              (void) sc_list_prepend (refine_list, elements[ci]);
            }
            /* We now recursively check the newly created elements for refinement. */
            t8_forest_adapt_refine_recursive (forest, ltree_id, tree->eclass, el_considered, scheme, refine_list,
                                              telements, &el_inserted, elements, &element_removed);
            el_coarsen = el_inserted;
          }
          else {
            (void) t8_element_array_push_count (telements, num_children);
            for (zz = 0; zz < num_children; zz++) {
              /* TODO: In a future version elements_from[zz] should be const and we should call t8_element_array_index_locidx (the const version). */
              elements[zz] = t8_element_array_index_locidx_mutable (telements, el_inserted + zz);
            }
            scheme->element_get_children (tree->eclass, elements_from[0], num_children, elements);
            el_inserted += (t8_locidx_t) num_children;
          }
          el_considered++;
        }
        else if (refine == -1) {
          /* The elements form a family and are to be coarsened. */
          /* Make room for one more new element. */
          elements[0] = t8_element_array_push (telements);
          /* Compute the parent of the current family.
           * This parent is now inserted in telements. */
          T8_ASSERT (scheme->element_get_level (tree->eclass, elements_from[0]) > 0);
          scheme->element_get_parent (tree->eclass, elements_from[0], elements[0]);
          /* num_siblings is now equivalent to the number of children of elements[0],
           * as num_siblings is always associated with elements_from*/
          num_children = num_siblings;
          el_inserted++;
          if (forest->set_adapt_recursive) {
            /* Adaptation is recursive.
             * We check whether the just generated parent is the last in its
             * family (and not the only one).
             * If so, we check this family for recursive coarsening. */
            const int child_id = scheme->element_get_child_id (tree->eclass, elements[0]);
            if (child_id > 0 && child_id == num_children - 1) {
              t8_forest_adapt_coarsen_recursive (forest, ltree_id, tree->eclass, el_considered, scheme, telements,
                                                 el_coarsen, &el_inserted, elements);
            }
          }
          el_considered += (t8_locidx_t) num_elements_to_adapt_callback;
        }
        else if (refine == 0) {
          /* The considered elements are neither to be coarsened nor is the first
           * one to be refined.
           * We copy the element to the new element array. */
          elements[0] = t8_element_array_push (telements);
          scheme->element_copy (tree->eclass, elements_from[0], elements[0]);
          el_inserted++;
          if (forest->set_adapt_recursive) {
            /* Adaptation is recursive.
             * If adaptation is recursive and this was the last element in its family
             * (and not the only one), we need to check for recursive coarsening. */
            const int child_id = scheme->element_get_child_id (tree->eclass, elements[0]);
            if (child_id > 0 && child_id == num_children - 1) {
              t8_forest_adapt_coarsen_recursive (forest, ltree_id, tree->eclass, el_considered, scheme, telements,
                                                 el_coarsen, &el_inserted, elements);
            }
          }
          el_considered++;
        }
        else {
          /* Remove the element */
          T8_ASSERT (refine == -2);
          element_removed = 1;
          el_considered++;
        }
      } /* End element loop */

      /* Check that if we had recursive adaptation, the refine list is now empty. */
      T8_ASSERT (!forest->set_adapt_recursive || refine_list->elem_count == 0);

      /* Set the new element offset of this tree */
      tree->elements_offset = el_offset;
      el_offset += el_inserted;
      /* Add to the new number of local elements. */
      forest->local_num_leaf_elements += el_inserted;
      /* Possibly shrink the telements array to the correct size */
      t8_element_array_resize (telements, el_inserted);

      /* It is not supported to delete all elements from a tree.
       * In this case, we will abort. */
      SC_CHECK_ABORTF (el_inserted != 0,
                       "ERROR: All elements of tree %i were removed. Removing all elements of a tree "
                       "is currently not supported. See also https://github.com/DLR-AMR/t8code/issues/1137.",
                       ltree_id);

      /* clean up */
      T8_FREE (elements);
      T8_FREE (elements_from);
    } /* End if (num_el_from > 0) */
  }   /* End tree loop */
  if (forest->set_adapt_recursive) {
    /* clean up */
    sc_list_destroy (refine_list);
  }

  /* We now adapted all local trees */
  /* Compute the new global number of elements */
  t8_forest_comm_global_num_leaf_elements (forest);

  /* Updating other processes about local (in)complete trees.
   * If the old forest already contained incomplete trees, 
   * this step is not necessary. */
  if (!forest_from->incomplete_trees) {
    T8_ASSERT (element_removed == 1 || element_removed == 0);
    int incomplete_trees;
    int mpiret = sc_MPI_Allreduce (&element_removed, &incomplete_trees, 1, sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
    SC_CHECK_MPI (mpiret);
    T8_ASSERT (incomplete_trees == 1 || incomplete_trees == 0);
    forest->incomplete_trees = incomplete_trees;
  }
  else {
    T8_ASSERT (forest_from->incomplete_trees == 1);
    forest->incomplete_trees = 1;
  }

  t8_global_productionf ("Done t8_forest_adapt with %lld total elements\n",
                         (long long) forest->global_num_leaf_elements);

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime += sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occurred on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End adapt %f %f\n", sc_MPI_Wtime (), forest->profile->adapt_runtime);
  }
}

T8_EXTERN_C_END ();

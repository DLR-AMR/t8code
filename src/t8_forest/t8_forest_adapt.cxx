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
#include <t8_forest.h>
#include <t8_data/t8_containers.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* Check the lastly inserted elements of an array for recursive coarsening.
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
t8_forest_adapt_coarsen_recursive (t8_forest_t forest, t8_locidx_t ltreeid,
                                   t8_locidx_t lelement_id,
                                   t8_eclass_scheme_c * ts,
                                   t8_element_array_t * telements,
                                   t8_locidx_t el_coarsen,
                                   t8_locidx_t * el_inserted,
                                   t8_element_t ** el_buffer)
{
  t8_element_t       *element, *parent;
  t8_element_t      **fam;
  t8_locidx_t         pos;
  size_t              elements_in_array;
  int                 num_children, i, isfamily;
  int                 child_id;
  /* el_inserted is the index of the last element in telements plus one.
   * el_coarsen is the index of the first element which could possibly
   * be coarsened. */

  elements_in_array = t8_element_array_get_count (telements);
  T8_ASSERT (*el_inserted == (t8_locidx_t) elements_in_array);
  T8_ASSERT (el_coarsen >= 0);
  element = t8_element_array_index_locidx (telements, *el_inserted - 1);
  /* TODO: This assumes that the number of children is the same for each
   *       element in that class. This may not be the case. */
  ts->t8_element_new (1, &parent);
  T8_ASSERT (ts->t8_element_level (element) > 0);
  ts->t8_element_parent (element, parent);
  num_children = ts->t8_element_num_children (parent);
  //T8_ASSERT (ts->t8_element_child_id (element) == num_children - 1);
  ts->t8_element_destroy (1, &parent);

  fam = el_buffer;
  pos = *el_inserted - num_children;
  isfamily = 1;
  child_id = ts->t8_element_child_id (element);
  while (isfamily && pos >= el_coarsen && child_id > 0 && child_id
         == num_children - 1) {
    isfamily = 1;
    /* Get all elements at indices pos, pos + 1, ... ,pos + num_children - 1 */
    for (i = 0; i < num_children; i++) {
      fam[i] = t8_element_array_index_locidx (telements, pos + i);
      if (ts->t8_element_child_id (fam[i]) != i) {
        /* These elements cannot form a family. Stop coarsening. */
        isfamily = 0;
        break;
      }
    }
    T8_ASSERT (!isfamily || ts->t8_element_is_family (fam));
    if (isfamily
        && forest->set_adapt_fn (forest, forest->set_from, ltreeid,
                                 lelement_id, ts, num_children, fam) < 0) {
      /* Coarsen the element */
      *el_inserted -= num_children - 1;
      /* remove num_children - 1 elements from the array */
      T8_ASSERT (elements_in_array == t8_element_array_get_count (telements));
      T8_ASSERT (ts->t8_element_level (element) > 0);
      ts->t8_element_parent (fam[0], fam[0]);
      num_children = ts->t8_element_num_children (fam[0]);
      elements_in_array -= num_children - 1;
      t8_element_array_resize (telements, elements_in_array);
      /* Set element to the new constructed parent. Since resizing the array
       * may change the position in memory, we have to do it after resizing. */
      element = t8_element_array_index_locidx (telements, pos);
    }
    else {
      /* If the elements are no family or
       * the family is not to be coarsened we abort the coarsening process */
      isfamily = 0;
    }
    pos -= num_children - 1;
  };
}

/* Check the lastly inserted element of an array for recursive refining.
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
                                  t8_eclass_scheme_c * ts,
                                  sc_list_t * elem_list,
                                  t8_element_array_t * telements,
                                  t8_locidx_t * num_inserted,
                                  t8_element_t ** el_buffer)
{
  t8_element_t       *insert_el;
  int                 num_children;
  int                 ci;

  if (elem_list->elem_count <= 0) {
    return;
  }

  while (elem_list->elem_count > 0) {
    /* Until the list is empty we
     * - remove the first element from the list.
     * - Check whether it should get refined
     * - If yes, we add all its children to the list
     * - If no, we add the element to the array of new elements
     */
    el_buffer[0] = (t8_element_t *) sc_list_pop (elem_list);
    num_children = ts->t8_element_num_children (el_buffer[0]);
    if (forest->set_adapt_fn (forest, forest->set_from, ltreeid, lelement_id,
                              ts, 1, el_buffer) > 0) {
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
    else {
      /* This element should not get refined,
       * we remove it from the buffer and add it to the array of new elements. */
      insert_el = t8_element_array_push (telements);
      ts->t8_element_copy (el_buffer[0], insert_el);
      ts->t8_element_destroy (1, el_buffer);
      (*num_inserted)++;
    }
  }
}

/* TODO: optimize this when we own forest_from */
void
t8_forest_adapt (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  sc_list_t          *refine_list = NULL;       /* This is only needed when we adapt recursively */
  t8_element_array_t *telements, *telements_from;
  t8_locidx_t         ltree_id, num_trees;
  t8_locidx_t         el_considered;
  t8_locidx_t         el_inserted;
  t8_locidx_t         el_coarsen;
  t8_locidx_t         num_el_from;
  t8_locidx_t         el_offset;
  size_t              num_children, zz, num_siblings,
    curr_num_siblings, curr_num_children;
  t8_tree_t           tree, tree_from;
  t8_eclass_scheme_c *tscheme;
  t8_element_t      **elements, **elements_from, *parent;
  int                 refine;
  int                 ci;
  unsigned int        subelement_type;
  unsigned int        num_subelements;
  long long int       count_subelements = 0;
#ifdef T8_ENABLE_DEBUG
  int                 is_family;
#endif

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
    num_children =
      tscheme->t8_element_num_children (t8_element_array_index_locidx
                                        (telements_from, 0));
    curr_num_children = num_children;
    curr_num_siblings =
      tscheme->t8_element_num_siblings (t8_element_array_index_locidx
                                        (telements_from, 0));
    tscheme->t8_element_new (1, &parent);
    /* Buffer for a family of new elements */
    elements = T8_ALLOC (t8_element_t *, num_children);
    /* Buffer for a family of old elements */
    elements_from = T8_ALLOC (t8_element_t *, curr_num_siblings);
    /* We now iterate over all elements in this tree and check them for refinement/coarsening/subelements. */
    while (el_considered < num_el_from) {
      int                 num_elements_to_adapt_fn;
#ifdef T8_ENABLE_DEBUG
      /* Will get set to 0 later if this is not a family */
      is_family = 1;
#endif

      /* Load the current element and at most num_siblings-1 many others into
       * the elements_from buffer. Stop when we are certain that they cannot from
       * a family.
       * At the end is_family will be true, if these elements form a family.
       */
      t8_element_t       *current_element =
        t8_element_array_index_locidx (telements_from, el_considered);

      num_siblings = tscheme->t8_element_num_siblings (current_element);

      if (num_siblings > curr_num_siblings) {
        elements_from =
          T8_REALLOC (elements_from, t8_element_t *, num_siblings);
        curr_num_siblings = num_siblings;
      }
      /* change: num_children into num_siblings */
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
        if ((size_t) tscheme->t8_element_child_id (elements_from[zz]) != zz) {
          break;
        }
      }

      if (zz != num_siblings
          || !tscheme->t8_element_is_family (elements_from)) {
        /* We are certain that the elements do not form a family.
         * So we will only pass the first element to the adapt callback. */
        num_elements_to_adapt_fn = 1;
#ifdef T8_ENABLE_DEBUG
        is_family = 0;
#endif
      }
      else {
        /* We will pass a family to the adapt callback */
        num_elements_to_adapt_fn = num_siblings;
      }
      T8_ASSERT (!is_family || tscheme->t8_element_is_family (elements_from));

      /* Pass the element, or the family to the adapt callback.
       * The output will be 
       *
       *      = -1 if we passed a family and it should get coarsened 
       *      =  0 if the element should remain as is 
       *      =  1 if the element should be refined (using the chosen recursive refinement scheme)
       *      >  1 if we use subelements
       * 
       * The values -1,0 and 1 will appear, if we use the standard refinement scheme of the given eclass.
       * 
       * Values above 1 will correspond to specific subelement types. 
       *  
       * For example the refine values for the 2D Quad scheme will be between -1 and 16. The values -1, 0 and 1 are for the standard refinement
       * and the values 2 to 16 correspond to the subelement types 1 to 15 (0001 to 1111 in base 2) and will be used by the element files of the quad 
       * scheme in order to remove hanging nodes.
       * 
       * It is up to the developer to use a reasonable range of subelement types for their use case. 
       * It is not allowed to use subelement types that intersect with already implemented ones.  
       * For the quad scheme, it would be reasonable to use subelement types >= 17, since 2 to 16 are already occupied by the hanging-
       * faces-removing ones. */

      refine = forest->set_adapt_fn (forest, forest->set_from, ltree_id,
                                     el_considered, tscheme,
                                     num_elements_to_adapt_fn, elements_from);

      /* We need to add a special condition for subelements that is independent of the refinement criteria:
       * After the adapt procedure, no subelements should be left in order to ensure that following modifications 
       * such as balance or remove_hanigng_faces work properly. 
       * Therefore, all subelements that "survive" the adaptation will be coarsened back to their parent quadrant. 
       * Note, that this is always valid for subelements in terms of the minimum level,
       * since subelements have the same level as their parent quadrant. */
      if (tscheme->t8_element_test_if_subelement (elements_from[0])
          && refine == 0) {
        refine = -1;
      }

#ifdef T8_ENABLE_DEBUG
      /* TODO: warning because of size_t instead of int */

      /* output for debugging */
      t8_debugf
        ("el_considered: %i/%i  refine: %i  is_family: %i  num_siblings: %i\n",
         el_considered + 1, num_el_from, refine, is_family, num_siblings);
#endif

      T8_ASSERT (is_family || refine >= 0);
      if (refine > 0
          && tscheme->t8_element_level (elements_from[0]) >=
          forest->maxlevel) {
        /* Only refine an element if it does not exceed the maximum level */
        refine = 0;
      }
      if (refine == 1) {
        /* The first element is to be refined using the standard quad scheme */
        num_children = tscheme->t8_element_num_children (elements_from[0]);
        if (num_children > curr_num_children) {
          elements = T8_REALLOC (elements, t8_element_t *, num_children);
          curr_num_children = num_children;
        }
        if (forest->set_adapt_recursive) {
          /* Create the children of this element */
          tscheme->t8_element_new (num_children, elements);     /*num_children or num_siblings? */
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
           * We can set this here to the next element after the current family, since a family that emerges from a refinement will never be coarsened */
          el_coarsen = el_inserted + num_children;
        }
        else {
          (void) t8_element_array_push_count (telements, num_children); /*num_children or num_siblings? */
          for (zz = 0; zz < num_children; zz++) {
            elements[zz] =
              t8_element_array_index_locidx (telements, el_inserted + zz);

          }
          tscheme->t8_element_children (elements_from[0], num_children,
                                        elements);
          el_inserted += num_children;
        }
        /* In case of a subelement, the parent quadrant is refined. 
         * Therfore, we can skip all subelement siblings as they are not needed anymore. */
        if (tscheme->t8_element_test_if_subelement (current_element)) {
          el_considered += num_siblings;
        }
        /* In case of a non-subelement element, we directly refine the quadrant and move on to the next element */
        else {
          el_considered++;
        }
      }
      else if (refine > 1) {    /* use subelements in this case */
        /* The subelement-callback function returns refine = subelement_type + 1 to avoid subelement_type = 1.
         * We can now undo this to use the subelement_type-values that match the binary encoding of the neighbour structure
         * (0001 should correspond to 1 and not to 2). */
        subelement_type = refine - 1;

        /* determing the number of subelements of the given type for memory allocation */
        num_subelements =
          tscheme->t8_element_get_number_of_subelements (subelement_type,
                                                         elements_from[0]);
        if (num_subelements > curr_num_children) {
          elements = T8_REALLOC (elements, t8_element_t *, num_subelements);
          curr_num_children = num_subelements;
        }

        (void) t8_element_array_push_count (telements, num_subelements);
        for (zz = 0; zz < num_subelements; zz++) {
          elements[zz] =
            t8_element_array_index_locidx (telements, el_inserted + zz);
        }
        tscheme->t8_element_to_subelement (elements_from[0], subelement_type,
                                           elements);
        el_inserted += num_subelements;
        el_considered++;

        /* Each time we entry this case, a parent element is refined into subelements.
         * We will count the global number of constructed subelements and give this number as additional output. */
        count_subelements += num_subelements;
      }
      else if (refine == -1) {

        /* The elements form a family and are to be coarsened. */
        /* Make room for one more new element. */
        elements[0] = t8_element_array_push (telements);
        /* Compute the parent of the current family.
         * This parent is now inserted in telements. */
        T8_ASSERT (tscheme->t8_element_level (elements_from[0]) > 0);
        tscheme->t8_element_parent (elements_from[0], elements[0]);
        el_inserted++;

        num_siblings = tscheme->t8_element_num_siblings (current_element);

        if (forest->set_adapt_recursive) {
          /* Adaptation is recursive.
           * We check whether the just generated parent is the last in its
           * family (and not the only one).
           * If so, we check this family for recursive coarsening. */
          const int           child_id = tscheme->t8_element_child_id (elements[0]);    /* elements[0] is the just constructed parent */
          if (child_id > 0 && (size_t) child_id == num_children - 1) {
            t8_forest_adapt_coarsen_recursive (forest, ltree_id,
                                               el_considered, tscheme,
                                               telements, el_coarsen,
                                               &el_inserted, elements);
          }
        }
        /* We coarsen a family of num_siblings many elements, skip them and move to the next element that is not part of this family */
        el_considered += num_siblings;
      }
      else {

        /* The considered elements are neither to be coarsened nor is the first
         * one to be refined.
         * We copy the element to the new element array. */
        T8_ASSERT (refine == 0);
        elements[0] = t8_element_array_push (telements);
        /* Note that it should not be allowed to copy a subelement at this point. 
         * Subelements should always be coarsened to their parent quadrant or their parent quadrant 
         * should be refined such that no subelements are left after calling the adapt function. */
        T8_ASSERT (!tscheme->t8_element_test_if_subelement (elements[0]));
        tscheme->t8_element_copy (elements_from[0], elements[0]);
        el_inserted++;
        const int           child_id =
          tscheme->t8_element_child_id (elements[0]);
        if (forest->set_adapt_recursive && child_id > 0
            && (size_t) tscheme->t8_element_child_id (elements[0])
            == num_children - 1) {
          /* If adaptation is recursive and this was the last element in its
           * family (and not the only one), we need to check for recursive coarsening. */
          t8_forest_adapt_coarsen_recursive (forest, ltree_id, el_considered,
                                             tscheme, telements, el_coarsen,
                                             &el_inserted, elements);
        }
        /* In this case the current element stays unchanged (gets copied) and we move on to the next element */
        el_considered++;
      }
    }                           /* end of while loop over elements */

    /* Check that if we had recursive adaptation, the refine list is now empty. */
    T8_ASSERT (!forest->set_adapt_recursive || refine_list->elem_count == 0);
    /* Set the new element offset of this tree */
    tree->elements_offset = el_offset;
    el_offset += el_inserted;
    /* Add to the new number of local elements. */
    forest->local_num_elements += el_inserted;
    /* Possibly shrink the telements array to the correct size */
    t8_element_array_resize (telements, el_inserted);

    /* clean up */
    T8_FREE (elements);
    T8_FREE (elements_from);
    tscheme->t8_element_destroy (1, &parent);
  }                             /* end of for loop over trees */
  if (forest->set_adapt_recursive) {
    /* clean up */
    sc_list_destroy (refine_list);
  }

  /* We now adapted all local trees */
  /* Compute the new global number of elements */
  t8_forest_comm_global_num_elements (forest);
  /* If any subelement is constructed, give output this number as an additional information. */
  if (count_subelements > 0) {
    t8_global_productionf
      ("Done t8_forest_adapt with %lld total elements, %lld of which are subelements.\n",
       (long long) forest->global_num_elements, count_subelements);
  }
  else {
    t8_global_productionf ("Done t8_forest_adapt with %lld total elements.\n",
                           (long long) forest->global_num_elements);
  }

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

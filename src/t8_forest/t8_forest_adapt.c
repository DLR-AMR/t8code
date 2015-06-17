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

/* TODO: optimize this when we own forest_from */
void
t8_forest_adapt (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  sc_list_t          *element_list;     /* This is only needed when we adapt recursive */
  sc_array_t         *telements, *telements_from;
  size_t              tt;
  t8_locidx_t         considered_elements;
  t8_locidx_t         inserted_elements;
  t8_topidx_t         treeid;
  size_t              num_children, zz;
  t8_tree_t           tree, tree_from;
  t8_eclass_scheme_t *tscheme;
  t8_element_t      **elements, **elements_from;
  int                 is_family;
  int                 refine;
  int                 level;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->set_adapt_recursive != -1);
  T8_ASSERT (forest->from_method == T8_FOREST_FROM_ADAPT);

  forest_from = forest->set_from;
  /* TODO: Allocate memory for the trees of forest.
   * Will we do this here or in an extra function? */
  T8_ASSERT (forest->trees->elem_count == forest_from->trees->elem_count);

  if (forest->set_adapt_recursive) {
    element_list = sc_list_new (NULL);
  }
  for (tt = 0; tt < forest->trees->elem_count; tt++) {
    treeid = tt + forest->first_local_tree;
    tree = (t8_tree_t) t8_sc_array_index_topidx (forest->trees, tt);
    tree_from = (t8_tree_t) t8_sc_array_index_topidx (forest_from->trees, tt);
    telements = &tree->elements;
    telements_from = &tree_from->elements;
    tscheme = forest->scheme->eclass_schemes[tree->eclass];
    considered_elements = 0;
    inserted_elements = 0;
    level = -1;
    /* TODO: this will generate problems with pyramidal elements */
    num_children = t8_eclass_num_children[tree->eclass];
    elements = T8_ALLOC (t8_element_t *, num_children);
    elements_from = T8_ALLOC (t8_element_t *, num_children);
    while (considered_elements < (t8_locidx_t) telements_from->elem_count) {
      is_family = 1;
      for (zz = 0; zz < num_children; zz++) {
        elements_from[zz] = t8_element_array_index (tscheme, telements_from,
                                                    considered_elements + zz);
        if ((size_t) t8_element_child_id (tscheme, elements_from[zz]) != zz) {
          elements_from[1] = NULL;
          is_family = 0;
          break;
        }
      }
      T8_ASSERT (!is_family || t8_element_is_family (tscheme, elements_from));
      refine = forest->set_adapt_fn (forest, treeid, tscheme, elements_from);
      if (refine > 0) {
        /* The first element is to be refined */
        if (forest->set_adapt_recursive) {
          /* TODO: recursive refinement */
        }
        else {
          /* add the children to the element array of the current tree */
          for (zz = 0; zz < num_children; zz++) {
            elements[zz] = (t8_element_t *) sc_array_push (telements);
            t8_element_child (tscheme, elements_from[0], zz, elements[zz]);
          }
          inserted_elements += num_children;
          considered_elements++;
        }
      }
      else if (refine < 0) {
        /* The elements form a family and are to be coarsened */
        if (forest->set_adapt_recursive) {
          /* TODO: recursive refinement */
        }
        else {
          elements[0] = (t8_element_t *) sc_array_push (telements);
          t8_element_parent (tscheme, elements_from[0], elements[0]);
          inserted_elements++;
          considered_elements += num_children;
        }
      }
      else {
        /* The considered elements are neither to be coarsened nor is the first
         * one to be refined */
        T8_ASSERT (refine == 0);
        if (forest->set_adapt_recursive) {
          if (t8_element_level (tscheme, elements_from[0]) > level + 1) {
            /* TODO: recursive refinement */
          }
          else if ((size_t) t8_element_child_id (tscheme, elements_from[0])
                   == num_children - 1) {
            /* TODO: recursiv refinement */
          }
        }
        else {
          elements[0] = (t8_element_t *) sc_array_push (telements);
          t8_element_copy (tscheme, elements_from[0], elements[0]);
          inserted_elements++;
          considered_elements++;
        }
      }
    }
    /* TODO: compute tree->element_offset */
  }
}

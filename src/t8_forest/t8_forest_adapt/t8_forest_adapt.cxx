/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/**
 * \file t8_forest_adapt/t8_forest_adapt.cxx
 * Implementation of the adaptation routine to refine and coarsen a forest of trees.
 */

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>


 void
 t8_forest_adapt_namespace::basic_adaptation::adapt()
 {
    T8_ASSERT (forest != nullptr);

    if (profiling) {
      profile_adaptation();
    }
    T8_ASSERT (forest_from != nullptr);

    collect_adapt_actions();

    /* Offset per tree in the source forest */ 
    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    /* Get the scheme used by the forest */
    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
        /* get the trees from both forests. */
        t8_tree_t tree = t8_forest_get_tree (forest, ltree_id);
        const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
        /* get the leaf arrays from both forests */
        t8_element_array_t *elements = &tree->leaf_elements;
        const t8_element_array_t *tree_elements_from = &tree_from->leaf_elements;
        /* Get the number of elements in the source tree */
        const t8_locidx_t num_el_from = (t8_locidx_t) t8_element_array_get_count (tree_elements_from);
        T8_ASSERT (num_el_from == t8_forest_get_tree_num_leaf_elements (forest_from, ltree_id));
        const t8_eclass_t tree_class = tree_from->eclass;
        /* Continue only if tree_from is not empty */
        if (num_el_from < 0){
            const t8_element_t *first_element_from = t8_element_array_index_locidx (tree_elements_from, 0);
            t8_locidx_t curr_size_elements_from = scheme->element_get_num_siblings (tree_class, first_element_from);
            /* index of the elements in source tree */
            t8_locidx_t el_considered = 0;
            /* index of the elements in target tree */
            t8_locidx_t el_inserted = 0;
            std::vector <const t8_element_t *> elements_temp;

            while (el_considered < num_el_from) {
                const t8_locidx_t num_siblings = scheme->element_get_num_siblings (tree_class, t8_element_array_index_locidx (tree_elements_from, el_considered));
                if (num_siblings > curr_size_elements_from) {
                    elements_temp.resize (num_siblings);
                    curr_size_elements_from = num_siblings;
                }
                for (int isibling = 0; isibling < num_siblings && el_considered + isibling < num_el_from; isibling++) {
                  elements_temp[isibling] = (const t8_element_t *) t8_element_array_index_locidx (tree_elements_from, el_considered + (t8_locidx_t )isibling);
                  if (scheme->element_get_child_id (tree_class, elements_temp[isibling]) != isibling) {
                    break;
                  }
                }
                const bool is_family = family_check (tree_elements_from, elements_temp, el_considered, scheme, tree_class);
                adapt_action action = adapt_actions[el_offset + el_considered];

                if (!is_family && action == COARSEN) {
                  action = KEEP;
                }
                /* Check that all siblings want to be coarsened */
                if (is_family && action == COARSEN) {
                  const auto start = adapt_actions.begin() + static_cast<size_t>(el_offset + el_considered);
                  const auto end = start + static_cast<size_t>(num_siblings);
                  if (!std::all_of(start, end, [](const adapt_action &a){ return a == COARSEN; })) {
                  action = KEEP;
                  }
                }

                switch (action) {
                case COARSEN:
                  el_inserted += manipulate_elements<COARSEN> (elements, tree_elements_from, scheme, tree_class,
                                                               el_inserted, el_offset + el_considered);
                  break;
                case KEEP:
                el_inserted += manipulate_elements<KEEP> (elements, tree_elements_from, scheme, tree_class,
                  el_inserted, el_offset + el_considered);
                  break;
                case REFINE:
                el_inserted += manipulate_elements<REFINE> (elements, tree_elements_from, scheme, tree_class,
                  el_inserted, el_offset + el_considered);
                  break;
                default:
                    {
                      t8_errorf ("Unknown adapt action.\n");
                      SC_ABORT_NOT_REACHED ();
                      break;
                    }
                }
                el_considered++;
            }
        }
        tree->elements_offset = el_offset;
        el_offset += num_el_from;

        
      }

 }
/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

#ifndef T8_FOREST_ADAPT_BATCHED_HXX
#define T8_FOREST_ADAPT_BATCHED_HXX

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>

/**
 * Collect adaptation actions for all elements in the source forest.
 */
struct batched_adapt_collector
{
  /**
   * Collect adaptation actions for all elements in the source forest.
   * \param [in] forest_from      The source forest from which to collect adaptation actions.
   * \param [out] adapt_actions   The vector to store the collected adaptation actions.
   * \param [in] callback         The callback function to determine the adaptation action for each element.
   * 
   * \warning Currently only used for testing purposes.
   */
  void
  collect_adapt_actions (const t8_forest_t forest_from, std::vector<t8_adapt::action> &adapt_actions,
                         batched_element_callback callback)
  {
    T8_ASSERT (forest_from != nullptr);

    t8_locidx_t el_offset = 0;
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest_from);
    const t8_locidx_t local_num_elements = t8_forest_get_local_num_leaf_elements (forest_from);
    adapt_actions.resize (local_num_elements);

    const t8_scheme *scheme = t8_forest_get_scheme (forest_from);

    /* For each element get the adaptation action */
    for (t8_locidx_t ltree_id = 0; ltree_id < num_trees; ltree_id++) {
      const t8_tree_t tree_from = t8_forest_get_tree (forest_from, ltree_id);
      const t8_eclass_t tree_class = tree_from->eclass;
      const t8_element_array_t *elements_from = &tree_from->leaf_elements;
      callback (forest_from, ltree_id, elements_from, scheme, tree_class, adapt_actions);
      el_offset += (t8_locidx_t) t8_element_array_get_count (elements_from);
    }
  }
};

#endif /* T8_FOREST_ADAPT_BATCHED_HXX */

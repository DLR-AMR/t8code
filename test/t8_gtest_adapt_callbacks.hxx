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

/** \file t8_gtest_adapt_callbacks.hxx
* Provide forest adapt callback functions that we use in our tests.
*/

#ifndef T8_GTEST_ADAPT_CALLBACKS
#define T8_GTEST_ADAPT_CALLBACKS

#include <src/t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx>

/** Adapt a forest such that always the first child of a
 * family is refined and no other elements. This results in a highly
 * imbalanced forest.
 *
 * This adapt callbacks requires an integer as forest user data.
 * This integer is the maximum refinement level.
 *
 * \param [in] forest       The forest to which the new elements belong.
 * \param [in] forest_from  The forest that is adapted.
 * \param [in] which_tree   The local tree containing \a elements.
 * \param [in] eclass       The eclass of \a which_tree.
 * \param [in] lelement_id  The local element id in \a forest_from in the tree of the current element.
 * \param [in] scheme       The scheme of the forest.
 * \param [in] is_family    If 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
 * \param [in] num_elements The number of entries in \a elements that are defined
 * \param [in] elements     Pointers to a family or, if \a is_family is zero,
 *                          pointer to one element.
 */
int
t8_test_adapt_first_child (const t8_forest_t forest, [[maybe_unused]]const t8_locidx_t which_tree, [[maybe_unused]]const t8_locidx_t lelement_id,
                         const t8_element_t *element, const t8_scheme_c *scheme, const t8_eclass_t tree_class);

/**
 * Adapt callback that refines every n-th global element, where \a n is given as template parameter.
 *
 * Optionally, an \a offset < \a n may be defined to not start with the first leaf element.
 *
 * \tparam n      every n-th global leaf element will be refined
 * \tparam offset global ID of the first element to refine
 *
 * Note: The argument list has to be the same as for \ref t8_forest_adapt_t, even
 *       if most arguments are unused. For their meaning, please refer to \ref t8_forest_adapt_t.
 *       (The following doxygen documentation is just to make sure it is technically documented.)
 *
 * \param[in] forest        "forest" argument of \ref t8_forest_adapt_t
 * \param[in] forest_from   "forest_from" argument of \ref t8_forest_adapt_t
 * \param[in] which_tree    "which_tree" argument of \ref t8_forest_adapt_t
 * \param[in] tree_class    "tree_class" argument of \ref t8_forest_adapt_t
 * \param[in] lelement_id   "lelement_id" argument of \ref t8_forest_adapt_t
 * \param[in] scheme        "scheme" argument of \ref t8_forest_adapt_t
 * \param[in] is_family     "is_family" argument of \ref t8_forest_adapt_t
 * \param[in] num_elements  "num_elements" argument of \ref t8_forest_adapt_t
 * \param[in] elements      "elements" argument of \ref t8_forest_adapt_t
 *
 * \return 1 if the element will be refined, 0 otherwise.
*/
template <int n, int offset = 0>
  requires (n > offset)
int
refine_every_nth_element_callback (const t8_forest_t forest,
                                   const t8_locidx_t which_tree,
                                   const t8_locidx_t lelement_id,
                                   [[maybe_unused]]const t8_element_t *element,
                                   [[maybe_unused]]const t8_scheme_c *scheme,
                                   [[maybe_unused]]const t8_eclass_t tree_class)
{
  // Compute global element ID (forest offset + tree offset + id in tree).
  const t8_forest_t forest_from = t8_forest_get_forest_from (forest);
  const t8_gloidx_t global_elem_id = t8_forest_get_first_local_leaf_element_id (forest_from)
                               + t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

  // Refine every n-th element, starting from offset (default: zero).
  return ((global_elem_id % n == offset) ? 1 : 0);
}

#endif /* T8_GTEST_ADAPT_CALLBACKS */

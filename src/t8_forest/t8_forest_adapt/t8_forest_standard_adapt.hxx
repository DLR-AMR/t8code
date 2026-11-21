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

/** \file t8_forest_standard_adapt.hxx
 * The adaptors for forests without subelements, with and without holes.
 */

#ifndef T8_FOREST_STANDARD_ADAPT_HXX
#define T8_FOREST_STANDARD_ADAPT_HXX

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>
#include <t8_forest/t8_forest_private.h>

/** Adaptor for a source forest without holes and without subelements. It uses all hooks of \ref t8_adaptor.
 * \tparam TRecursive If true, adapt recursively, see \ref t8_adaptor.
 */
template <bool TRecursive>
struct t8_standard_adaptor: t8_adaptor<t8_standard_adaptor<TRecursive>, TRecursive>
{
  using t8_adaptor<t8_standard_adaptor, TRecursive>::t8_adaptor;
};

/** Adaptor for a source forest with holes (incomplete trees) and without subelements.
 * Families may be incomplete, so they are detected with \ref t8_forest_is_incomplete_family and the callback
 * receives incomplete families as well.
 * \tparam TRecursive If true, adapt recursively, see \ref t8_adaptor.
 */
template <bool TRecursive>
struct t8_incomplete_adaptor: t8_adaptor<t8_incomplete_adaptor<TRecursive>, TRecursive>
{
  using t8_adaptor<t8_incomplete_adaptor, TRecursive>::t8_adaptor;

  /** Collect the (possibly incomplete) family starting at \a context.el_considered.
   * \see t8_adaptor::get_family
   * \param [in] context The current tree.
   * \return The family starting at \a context.el_considered, or this element alone.
   */
  t8_adapt_family
  get_family (const t8_adapt_tree_context &context)
  {
    t8_adapt_family family;
    const int num_siblings = context.scheme->element_get_num_siblings (
      context.tree_class, t8_element_array_index_locidx (context.elements_from, context.el_considered));
    T8_ASSERT (num_siblings <= T8_ECLASS_MAX_CHILDREN);
    int num_collected = 0;
    for (; num_collected < num_siblings && context.el_considered + num_collected < context.num_elements_from;
         num_collected++) {
      family.elements[num_collected]
        = t8_element_array_index_locidx_mutable (context.elements_from, context.el_considered + num_collected);
    }
    const int family_size = t8_forest_is_incomplete_family (this->forest_from, context.ltree_id, context.el_considered,
                                                            family.elements.data (), num_collected);
    family.is_family = family_size > 0;
    family.num_elements = family.is_family ? family_size : 1;
    return family;
  }

  /** Collect the elements of the new tree from \a first to its end as an (possibly incomplete) family.
   * \see t8_adaptor::get_new_family
   * \param [in] context The current tree.
   * \param [in] first   The index of the first family member in \a context.elements, as found by
   *                     \ref t8_adaptor::find_family_start.
   * \return The family from \a first to the end of \a context.elements.
   */
  t8_adapt_family
  get_new_family (const t8_adapt_tree_context &context, const t8_locidx_t first)
  {
    t8_adapt_family family;
    family.num_elements = (int) t8_element_array_get_count (context.elements) - first;
    T8_ASSERT (0 < family.num_elements && family.num_elements <= T8_ECLASS_MAX_CHILDREN);
    for (int ielement = 0; ielement < family.num_elements; ielement++) {
      family.elements[ielement] = t8_element_array_index_locidx_mutable (context.elements, first + ielement);
    }
    family.is_family = true;
    return family;
  }
};

#endif /* T8_FOREST_STANDARD_ADAPT_HXX */

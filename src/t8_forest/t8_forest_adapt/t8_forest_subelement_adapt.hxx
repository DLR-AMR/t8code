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

/** \file t8_forest_subelement_adapt.hxx
 * The adaptor for complete forests with a subelement scheme.
 */

#ifndef T8_FOREST_SUBELEMENT_ADAPT_HXX
#define T8_FOREST_SUBELEMENT_ADAPT_HXX

#include <t8_forest/t8_forest_adapt/t8_forest_adapt.hxx>
#include <t8_schemes/t8_subelement/t8_subelement.hxx>

/** Adaptor for a complete source forest with a subelement scheme.
 * An adapt callback return value greater than one requests the transition of the element into the subelements of
 * type return value minus one. The offset avoids a collision of subelement type one with regular refinement.
 * \tparam TRecursive If true, adapt recursively, see \ref t8_adaptor. Subelements are only created for elements
 *                    of the source forest, not for elements created by recursive refinement.
 */
template <bool TRecursive>
struct t8_subelement_adaptor: t8_adaptor<t8_subelement_adaptor<TRecursive>, TRecursive>
{
  using t8_adaptor<t8_subelement_adaptor, TRecursive>::t8_adaptor;

  /** Replace the first element of \a family by its subelements of type \a type.
   * \see t8_adaptor::refine_special
   * \param [in, out] context The current tree.
   * \param [in]      family  The elements passed to the callback.
   * \param [in]      type    The subelement type.
   * \return The number of consumed source elements.
   */
  t8_locidx_t
  refine_special (t8_adapt_tree_context &context, const t8_adapt_family &family, const int type)
  {
    T8_ASSERT (t8_eclass_scheme_is_subelement (context.scheme, context.tree_class));
    const int num_subelements = context.scheme->element_get_num_children (context.tree_class, family.elements[0], type);
    (void) t8_element_array_push_count (context.elements, num_subelements);
    context.scheme->element_get_children (context.tree_class, family.elements[0], num_subelements,
                                          this->last_elements (context, num_subelements).data (), type);
    return 1;
  }
};

#endif /* T8_FOREST_SUBELEMENT_ADAPT_HXX */

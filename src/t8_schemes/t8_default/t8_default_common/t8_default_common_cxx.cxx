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

#include <sc_functions.h>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/** Compute the number of corners of a given element. */
int
t8_default_scheme_common_c::t8_element_num_corners (const t8_element_t *elem) const
{
  /* use the lookup table of the eclasses.
   * Pyramids should implement their own version of this function. */
  return t8_eclass_num_vertices[eclass];
}

t8_element_shape_t
t8_default_scheme_common_c::t8_element_shape (const t8_element_t *elem) const
{
  return eclass;
}

/* Given an element's level and dimension, return the number of leaves it
 * produces at a given uniform refinement level */
static inline t8_gloidx_t
count_leaves_from_level (int element_level, int refinement_level, int dimension)
{
  return element_level > refinement_level ? 0 : sc_intpow64 (2, dimension * (refinement_level - element_level));
}

t8_gloidx_t
t8_default_scheme_common_c::t8_element_count_leaves (const t8_element_t *t, int level) const
{

  int element_level = t8_element_level (t);
  t8_element_shape_t element_shape;
  int dim = t8_eclass_to_dimension[eclass];
  element_shape = t8_element_shape (t);
  if (element_shape == T8_ECLASS_PYRAMID) {
    int level_diff = level - element_level;
    return element_level > level ? 0 : 2 * sc_intpow64 (8, level_diff) - sc_intpow64 (6, level_diff);
  }
  return count_leaves_from_level (element_level, level, dim);
}

/* Count the number of siblings.
 * The number of children is 2^dim for each element, except for pyramids.
 * TODO: For pyramids we will have to implement a standalone version in the pyramid scheme. */
int
t8_default_scheme_common_c::t8_element_num_siblings (const t8_element_t *elem) const
{
  const int dim = t8_eclass_to_dimension[eclass];
  T8_ASSERT (eclass != T8_ECLASS_PYRAMID);
  return sc_intpow (2, dim);
}

t8_gloidx_t
t8_default_scheme_common_c::t8_element_count_leaves_from_root (int level) const
{
  if (eclass == T8_ECLASS_PYRAMID) {
    return 2 * sc_intpow64u (8, level) - sc_intpow64u (6, level);
  }
  int dim = t8_eclass_to_dimension[eclass];
  return count_leaves_from_level (0, level, dim);
}

T8_EXTERN_C_END ();

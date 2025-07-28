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

#ifndef T8_STANDALONE_ELEMENTS_HXX
#define T8_STANDALONE_ELEMENTS_HXX

#include <t8.h>
#include <array>
#include <bitset>

#define t8_standalone_element t8_standalone

/** Dimension of the standalone element types */
constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 3 };

/** Maximum level of the standalone element types
 * \note The maxlevel is lower than 255 so that we can use \ref t8_element_level (uint8_t)
 * to iterate to maxlevel:
 * for (t8_element_level level = 0; level <= T8_ELEMENT_MAXLEVEL[T8_ECLASS_VERTEX]; ++level)
 * Otherwise, t8_element_level would overflow after 255 and we would have an infinite loop.
*/
constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = { 254, 30, 30, 29, 21, 21, 21, 18 };

/** Maximum number of faces of the standalone element types */
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_COUNT] = { 1, 2, 4, 3, 6, 4, 5, 5 };

/** Number of children of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_COUNT] = { 1, 2, 4, 4, 8, 8, 8, 10 };

/** Number of corners (vertices) of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_CORNERS[T8_ECLASS_COUNT] = { 1, 2, 4, 3, 8, 4, 6, 5 };

/** Actual number of faces of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_FACES[T8_ECLASS_COUNT] = { 0, 2, 4, 3, 6, 4, 5, 5 };

/** Number of facechildren of the standalone element types */
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACECHILDREN[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 4, 4, 4, 4 };

constexpr uint8_t T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_COUNT] = { 0, 0, 0, 1, 0, 3, 1, 2 };

typedef int32_t t8_element_coord;
typedef uint8_t t8_element_level;
typedef uint8_t t8_cube_id;
typedef uint8_t t8_child_id;

template <t8_eclass_t TEclass>
using t8_element_type = std::bitset<T8_ELEMENT_NUM_EQUATIONS[TEclass]>;

template <t8_eclass_t TEclass>
using t8_element_coords = std::array<t8_element_coord, T8_ELEMENT_DIM[TEclass]>;
template <t8_eclass_t TEclass>
struct t8_standalone_element
{
  /** The coordinates of the anchor vertex of the element. */
  t8_element_coords<TEclass> coords;
  /** The refinement level of the element relative to the root at level 0. */
  t8_element_level level;
};

#endif /* T8_STANDALONE_ELEMENTS_HXX */

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

constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 3 };
constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = { 255, 30, 30, 29, 21, 21, 21, 18 };
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_COUNT] = { 1, 2, 4, 3, 6, 4, 5, 5 };
constexpr uint8_t T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_COUNT] = { 1, 2, 4, 4, 8, 8, 8, 10 };
constexpr uint8_t T8_ELEMENT_NUM_CORNERS[T8_ECLASS_COUNT] = { 1, 2, 4, 3, 8, 4, 6, 5 };
constexpr uint8_t T8_ELEMENT_NUM_FACES[T8_ECLASS_COUNT] = { 0, 2, 4, 3, 6, 4, 5, 5 };
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
struct t8_standalone_element : t8_element_base<t8_standalone_scheme<TEclass>>
{
  /** The coordinates of the anchor vertex of the element. */
  t8_element_coords<TEclass> coords;
  /** The refinement level of the element relative to the root at level 0. */
  t8_element_level level;
};



#endif /* T8_STANDALONE_ELEMENTS_HXX */

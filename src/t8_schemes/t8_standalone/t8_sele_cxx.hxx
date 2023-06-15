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

#ifndef T8_SELE_HXX
#define T8_SELE_HXX

#include <t8.h>
#include <array>
#include <bitset>

// name is to long.. maybe t8_sele_t or t8_selement_t or t8_element_t 
#define t8_standalone_element_t t8_sele_t

typedef uint32_t t8_element_coord_t;
typedef uint8_t  t8_element_level_t;
typedef int8_t   t8_cube_id_t;

constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = {0, 1, 2, 2, 3, 3, 3, 3};
constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {255, 30, 29, 21, 18, 21, 21, 21};
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_COUNT] = {1, 2, 4, 3, 6, 4, 5, 5};
constexpr uint8_t T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_COUNT] = {0, 2, 4, 4, 8, 8, 8, 10};
constexpr uint8_t T8_ELEMENT_NUM_CORNERS[T8_ECLASS_COUNT] = {1, 2, 4, 3, 8, 4, 6, 5};
constexpr uint8_t T8_ELEMENT_NUM_FACES[T8_ECLASS_COUNT] = {0, 2, 4, 3, 6, 4, 5, 5};

/* To be filled */
constexpr uint8_t T8_ELEMENT_TYPE[T8_ECLASS_COUNT] = {0, 0, 0, 0, 0, 0, 0, 0};
constexpr uint8_t T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_COUNT] = {0, 0, 0, 0, 0, 0, 0, 0};

/** The type of pyramid in 0, ...,7. The first 6 types describe tetrahedra.
 * Type 6 is an upward facing pyramid.
 * Type 7 is a downward facing pyramid.
*/


template <t8_eclass_t eclass_T>
struct t8_standalone_element_t
{
  /** The refinement level of the element relative to the root at level 0. */
  t8_element_level_t level; 

  /** Bit array: which inequality is fulfilled at which level. */
  std::bitset<T8_ELEMENT_NUM_EQUATIONS[eclass_T]>  type;
  std::array<t8_element_coord_t, T8_ELEMENT_DIM[eclass_T]> coords;
};

#endif /* T8_SELE_HXX */

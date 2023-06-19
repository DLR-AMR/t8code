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

#ifndef T8_SELE_LUT_TRIANGLE_HXX
#define T8_SELE_LUT_TRIANGLE_HXX
// TRIANGLES

template<>
constexpr int8_t t8_element_type_Iloc_to_childtype<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TRIANGLE]] = {
  {0, 0, 1, 0},
  {1, 0, 1, 1}
};

template<>
constexpr int8_t t8_element_type_Iloc_to_childcubeid<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]]
  [T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TRIANGLE]] = {
  {0, 1, 1, 3},
  {0, 2, 2, 3}
};

template<>
constexpr int8_t t8_element_type_cubeid_to_parenttype<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][1<<T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]
  = {
  {0, 0, 1, 0},
  {1, 0, 1, 1}
};

template<>
constexpr int8_t t8_element_type_cubeid_to_Iloc<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][1<<T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]
  = {
  {0, 1, 1, 3},
  {0, 2, 2, 3}
};

template<>
constexpr int8_t t8_type_edge_equations<T8_ECLASS_TRIANGLE>[T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][2] = {
  {0, 1}
};

template<>
constexpr int8_t t8_type_vertex_dim_to_binary<T8_ECLASS_TRIANGLE>[1 << T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]]
  [T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]] = {
  {
   {0, 0},
   {1, 0},
   {1, 1}
   },
  {
  /* PAY ATTENTION TO ORDERING OF ELEMENTS! What is the ordering? **/
   {0, 0},
   {0, 1},
   {1, 1}
   }
};

#endif /* T8_SELE_LUT_TRIANGLE_HXX */
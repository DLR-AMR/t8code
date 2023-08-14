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

constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = {0, 1, 2, 2, 3, 3, 3, 3};
constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {255, 30, 29, 29, 21, 21, 21, 18};
//constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {255, 15, 15, 15, 15, 15, 15, 15};
//constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {3, 3, 3, 3, 3, 3, 3, 3};
//constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {5, 5, 5, 5, 5, 5, 5, 5};
//constexpr uint8_t T8_ELEMENT_MAXLEVEL[T8_ECLASS_COUNT] = {2, 2, 2, 2, 2, 2, 2, 2};
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_COUNT] = {1, 2, 4, 3, 6, 4, 5, 5};
constexpr uint8_t T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_COUNT] = {1, 2, 4, 4, 8, 8, 8, 10};
constexpr uint8_t T8_ELEMENT_NUM_CORNERS[T8_ECLASS_COUNT] = {1, 2, 4, 3, 8, 4, 6, 5};
constexpr uint8_t T8_ELEMENT_NUM_FACES[T8_ECLASS_COUNT] = {0, 2, 4, 3, 6, 4, 5, 5};
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACECHILDREN[T8_ECLASS_COUNT] = {0, 1, 2, 2, 4, 4, 4, 4};

constexpr uint8_t T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_COUNT] = {0, 0, 0, 1, 0, 3, 1, 2};

/**PARENT CHILD BIJECTION*/
template<t8_eclass_t eclass_T>
constexpr int8_t t8_element_type_Iloc_to_childtype[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_CHILDREN[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_element_type_Iloc_to_childcubeid[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_CHILDREN[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_element_type_cubeid_to_parenttype[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][1<<T8_ELEMENT_DIM[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_element_type_cubeid_to_Iloc[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][1<<T8_ELEMENT_DIM[eclass_T]];


/**TYPE EQUATIONS*/
template<t8_eclass_t eclass_T>
constexpr int8_t t8_type_edge_equations[T8_ELEMENT_NUM_EQUATIONS[eclass_T]][2];


/**VERTEX*/
template<t8_eclass_t eclass_T>
constexpr int8_t t8_type_vertex_dim_to_binary[1 << T8_ELEMENT_NUM_EQUATIONS[eclass_T]]
  [T8_ELEMENT_NUM_CORNERS[eclass_T]][T8_ELEMENT_DIM[eclass_T]];

/**FACE*/
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_face_internal[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_typebit[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_sign[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_facenormal_dim[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_neighface[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_cubeid_face_to_parentface[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][1<<T8_ELEMENT_DIM[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_is_1_boundary[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_last_facechilds_cubeid[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_facechildid_to_childid[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]][T8_ELEMENT_MAX_NUM_FACECHILDREN[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_childid_face_to_childface[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_CHILDREN[eclass_T]][T8_ELEMENT_MAX_NUM_FACECHILDREN[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_face_to_tree_face[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_rootface_dim_to_facedim[T8_ELEMENT_NUM_FACES[eclass_T]][T8_ELEMENT_DIM[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_rootface_eq_to_faceeq[T8_ELEMENT_NUM_FACES[eclass_T]][T8_ELEMENT_NUM_EQUATIONS[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr int8_t t8_sele_lut_type_rootface_to_face[1<<T8_ELEMENT_NUM_EQUATIONS[eclass_T]][T8_ELEMENT_NUM_FACES[eclass_T]];
template<t8_eclass_t eclass_T>
constexpr t8_eclass_t t8_sele_lut_rootface_to_eclass[T8_ELEMENT_NUM_FACES[eclass_T]];



#include "t8_sele_lut_triangle_cxx.hxx"
#include "t8_sele_lut_prism_cxx.hxx"
#include "t8_sele_lut_pyra_cxx.hxx"
#include "t8_sele_lut_tet_cxx.hxx"


typedef uint32_t t8_element_coord_t;
typedef uint8_t  t8_element_level_t;
typedef int8_t   t8_cube_id_t;

template <t8_eclass_t eclass_T>
using t8_element_type_t = std::bitset<T8_ELEMENT_NUM_EQUATIONS[eclass_T]>;

template <t8_eclass_t eclass_T>
using t8_element_coords_t = std::array<t8_element_coord_t, T8_ELEMENT_DIM[eclass_T]>;

template <t8_eclass_t eclass_T>
struct t8_standalone_element_t
{
  /** The refinement level of the element relative to the root at level 0. */
  t8_element_level_t level; 

  /** The coordinates of the anchor vertex of the element. */
  t8_element_coords_t<eclass_T> coords;

  /** Bit array: which inequality is fulfilled at which level. */
  t8_element_type_t<eclass_T> type;
};

#endif /* T8_SELE_HXX */

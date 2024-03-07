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

#ifndef T8_SELE_LUT_TET_HXX
#define T8_SELE_LUT_TET_HXX
// TETS

// clang-format off
template<>
constexpr int8_t t8_element_type_Iloc_to_childtype<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TET]] = {
  {0, 0, 1, 5, 0, 2, 6, 0},
  {1, 0, 1, 2, 1, 5, 7, 1},
  {2, 2, 6, 7, 0, 1, 2, 2},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {5, 5, 6, 7, 0, 1, 5, 5},
  {6, 0, 2, 6, 5, 6, 7, 6},
  {7, 1, 5, 7, 2, 6, 7, 7}
};

template<>
constexpr int8_t t8_element_type_Iloc_to_childcubeid<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]]
  [T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TET]] = {
  {0, 1, 1, 1, 3, 3, 3, 7},
  {0, 2, 2, 2, 3, 3, 3, 7},
  {0, 1, 1, 1, 5, 5, 5, 7},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {0, 2, 2, 2, 6, 6, 6, 7},
  {0, 4, 4, 4, 5, 5, 5, 7},
  {0, 4, 4, 4, 6, 6, 6, 7},
};

template<>
constexpr int8_t t8_element_type_cubeid_to_parenttype<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][1<<T8_ELEMENT_DIM[T8_ECLASS_TET]]
  = {
  {0, 0, 1, 0, 6, 2, 5, 0},
  {1, 0, 1, 1, 7, 2, 5, 1},
  {2, 2, 1, 0, 6, 2, 7, 2},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {5, 0, 5, 1, 7, 6, 5, 5},
  {6, 2, 5, 0, 6, 6, 7, 6},
  {7, 2, 5, 1, 7, 6, 7, 7}
};

template<>
constexpr int8_t t8_element_type_cubeid_to_Iloc<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][1<<T8_ELEMENT_DIM[T8_ECLASS_TET]]
  = {
  {0, 1, 1, 4, 1, 4, 4, 7},
  {0, 2, 2, 4, 1, 5, 5, 7},
  {0, 1, 3, 5, 2, 6, 4, 7},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1, -1, -1, -1},
  {0, 3, 1, 5, 2, 4, 6, 7},
  {0, 2, 2, 6, 3, 5, 5, 7},
  {0, 3, 3, 6, 3, 6, 6, 7}
};



/**corresponds to x >= y >= z*/
template<>
constexpr int8_t t8_type_edge_equations<T8_ECLASS_TET>[T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][2] = {
  {0, 1},
  {1, 2},
  {0, 2}
};

template<>
constexpr int8_t t8_type_vertex_dim_to_binary<T8_ECLASS_TET>[1 << T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]]
  [T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TET]][T8_ELEMENT_DIM[T8_ECLASS_TET]] = {
  {
   {0, 0, 0},
   {1, 0, 0},
   {1, 1, 0},
   {1, 1, 1}
   },
  {
   {0, 0, 0},
   {0, 1, 0},
   {1, 1, 0},
   {1, 1, 1}
   },
  {
   {0, 0, 0},
   {1, 0, 0},
   {1, 0, 1},
   {1, 1, 1}
   },
  {
   {-1, -1, -1},
   {-1, -1, -1},
   {-1, -1, -1},
   {-1, -1, -1}
   },
  {
   {-1, -1, -1},
   {-1, -1, -1},
   {-1, -1, -1},
   {-1, -1, -1}
   },
  {
   {0, 0, 0},
   {0, 1, 0},
   {0, 1, 1},
   {1, 1, 1}
   },
  {
   {0, 0, 0},
   {0, 0, 1},
   {1, 0, 1},
   {1, 1, 1}
   },
  {
   {0, 0, 0},
   {0, 0, 1},
   {0, 1, 1},
   {1, 1, 1}
   }
};
template<>
constexpr int8_t t8_sele_lut_face_internal<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]] = {
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0},
  {0,1,1,0}
};
template<>
constexpr int8_t t8_sele_lut_type_face_to_typebit<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]] = {
  {-1, 0, 1,-1},
  {-1, 0, 2,-1},
  {-1, 2, 1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1, 1, 2,-1},
  {-1, 2, 0,-1},
  {-1, 1, 0,-1}
};

template<>
constexpr int8_t t8_sele_lut_type_face_to_sign<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]] = {
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1},
  { 1, 0, 0,-1}
};

template<>
constexpr int8_t t8_sele_lut_type_face_to_facenormal_dim<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]] = {
  { 0,-1,-1, 2},
  { 1,-1,-1, 2},
  { 0,-1,-1, 1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  { 1,-1,-1, 0},
  { 2,-1,-1, 1},
  { 2,-1,-1, 0}
};

template<>
constexpr int8_t t8_sele_lut_type_face_to_is_1_boundary<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0},
  {1,0,0,0}
};


template<>
constexpr int8_t t8_sele_lut_type_face_to_neighface<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
  { 3, 1, 2, 0},
  { 3, 1, 2, 0},
  { 3, 1, 2, 0},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  { 3, 1, 2, 0},
  { 3, 1, 2, 0},
  { 3, 1, 2, 0}
};


/** TODO! Fill*/
template<>
constexpr int8_t t8_sele_lut_type_face_to_tree_face<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
  { 0, 1, 2, 3},
  {-1,-1,-1, 3},
  { 0,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1, 2,-1,-1},
  {-1,-1, 1,-1},
  {-1,-1,-1,-1}
};


template<>
constexpr int8_t t8_sele_lut_type_cubeid_face_to_parentface<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][1<<T8_ELEMENT_DIM[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
{
  {-1,1,2,3},
  {0,-1,2,3},
  {-1,-1,-1,3},
  {0,1,-1,3},
  {-1,2,-1,-1},
  {0,-1,-1,-1},
  {-1,-1,1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,-1,-1,3},
  {0,-1,2,3},
  {0,1,-1,3},
  {-1,2,-1,-1},
  {-1,-1,1,-1},
  {0,-1,-1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {0,-1,2,3},
  {-1,2,-1,-1},
  {0,-1,-1,-1},
  {-1,-1,-1,3},
  {0,1,-1,3},
  {-1,-1,1,-1},
  {0,1,2,-1}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
   {-1,1,2,3},
  {-1,2,-1,-1},
  {0,-1,2,3},
  {0,-1,-1,-1},
  {-1,-1,-1,3},
  {-1,-1,1,-1},
  {0,1,-1,3},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,-1,-1,3},
  {-1,2,-1,-1},
  {-1,-1,1,-1},
  {0,-1,2,3},
  {0,1,-1,3},
  {0,-1,-1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,2,-1,-1},
  {-1,-1,-1,3},
  {-1,-1,1,-1},
  {0,-1,2,3},
  {0,-1,-1,-1},
  {0,1,-1,3},
  {0,1,2,-1}
}
};

template<>
constexpr int8_t t8_sele_lut_type_childid_face_to_childface<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
{
  {-1,1,2,3},
  {0,-1,2,3},
  {-1,-1,-1,3},
  {-1,-1,1,-1},
  {0,1,-1,3},
  {0,-1,-1,-1},
  {-1,2,-1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,-1,-1,3},
  {0,-1,2,3},
  {-1,-1,1,-1},
  {0,1,-1,3},
  {0,-1,-1,-1},
  {-1,2,-1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {0,-1,2,3},
  {-1,-1,-1,3},
  {-1,-1,1,-1},
  {0,-1,-1,-1},
  {-1,2,-1,-1},
  {0,1,-1,3},
  {0,1,2,-1}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
  {-1,1,2,3},
  {0,-1,2,3},
  {-1,-1,1,-1},
  {-1,-1,-1,3},
  {-1,2,-1,-1},
  {0,-1,-1,-1},
  {0,1,-1,3},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,-1,1,-1},
  {-1,-1,-1,3},
  {0,-1,2,3},
  {-1,2,-1,-1},
  {0,1,-1,3},
  {0,-1,-1,-1},
  {0,1,2,-1}
},
{
  {-1,1,2,3},
  {-1,-1,1,-1},
  {-1,-1,-1,3},
  {0,-1,2,3},
  {-1,2,-1,-1},
  {0,-1,-1,-1},
  {0,1,-1,3},
  {0,1,2,-1}
}
};

template<>
constexpr int8_t t8_sele_lut_type_face_to_last_facechilds_cubeid<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]]={
  {7,7,7,3},
  {7,7,7,3},
  {7,7,7,5},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {7,7,7,6},
  {7,7,7,5},
  {7,7,7,6}
};

template<>
constexpr int8_t t8_sele_lut_type_face_facechildid_to_childid<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TET]][T8_ELEMENT_MAX_NUM_FACECHILDREN[T8_ECLASS_TET]]={
{
  {1,4,5,7},
  {0,4,6,7},
  {0,1,3,7},
  {0,1,2,4}
},
{
  {2,4,5,7},
  {0,4,6,7},
  {0,2,3,7},
  {0,1,2,4}
},
{
  {1,4,6,7},
  {0,5,6,7},
  {0,1,3,7},
  {0,1,2,6}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1}
},
{
  {1,5,6,7},
  {0,4,6,7},
  {0,1,2,7},
  {0,1,3,6}
},
{
  {3,5,6,7},
  {0,4,5,7},
  {0,1,3,7},
  {0,2,3,5}
},
{
  {3,5,6,7},
  {0,4,6,7},
  {0,1,3,7},
  {0,2,3,6}
}
};


template<>
constexpr int8_t t8_sele_lut_rootface_dim_to_facedim<T8_ECLASS_TET>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]][T8_ELEMENT_DIM[T8_ECLASS_TET]]={
  {-1,0,1},
  {0,0,1},
  {0,1,1},
  {0,1,-1}
};
template<>
constexpr int8_t t8_sele_lut_rootface_eq_to_faceeq<T8_ECLASS_TET>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]][T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]]={
  {-1, 0,-1},
  {-1, 0, 0},
  { 0,-1, 0},
  { 0,-1,-1}  
};

template<>
constexpr int8_t t8_sele_lut_type_rootface_to_face<T8_ECLASS_TET>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TET]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]]={
  { 0, 1, 2, 3},
  {-1,-1,-1, 3},
  { 0,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1,-1,-1},
  {-1,-1, 1,-1},
  {-1, 2,-1,-1},
  {-1,-1,-1,-1} 
};

template<>
constexpr t8_eclass_t t8_sele_lut_rootface_to_eclass<T8_ECLASS_TET>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TET]]={
  T8_ECLASS_TRIANGLE,
  T8_ECLASS_TRIANGLE,
  T8_ECLASS_TRIANGLE,
  T8_ECLASS_TRIANGLE
};
// clang-format on
#endif /* T8_SELE_LUT_TET_HXX */
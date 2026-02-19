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

#ifndef T8_SELE_LUT_TRIANGLE_HXX
#define T8_SELE_LUT_TRIANGLE_HXX
// TRIANGLES
// clang-format off
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


template<>
constexpr int8_t t8_standalone_lut_face_internal<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]] = {
  {0,1,0},
  {0,1,0}
};
template<>
constexpr int8_t t8_standalone_lut_type_face_to_typebit<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]] = {
  {-1,0,-1},
  {-1,0,-1}
};

template<>
constexpr int8_t t8_standalone_lut_type_face_to_sign<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]] = {
  {1,0,-1},
  {1,0,-1}
};

template<>
constexpr int8_t t8_standalone_lut_type_face_to_is_1_boundary<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  {1,0,0},
  {1,0,0}
};

template<>
constexpr int8_t t8_standalone_lut_type_face_to_facenormal_dim<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]] = {
  {0,-1,1},
  {1,-1,0}
};


template<>
constexpr int8_t t8_standalone_lut_type_face_to_neighface<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  {2,1,0},
  {2,1,0}
};

template<>
constexpr int8_t t8_standalone_lut_type_face_to_tree_face<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  {0,1,2},
  {-1,-1,-1}
};


template<>
constexpr int8_t t8_standalone_lut_type_cubeid_face_to_parentface<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][1<<T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
{
  {-1,1,2},
  {0,-1,2},
  {-1,-1,-1},
  {0,1,-1}
},
{
  {-1,1,2},
  {-1,-1,-1},
  {0,-1,2},
  {0,1,-1}
}
};

template<>
constexpr int8_t t8_standalone_lut_type_childid_face_to_childface<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_CHILDREN[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
{
  {-1,1,2},
  {0,-1,2},
  {-1,-1,-1},
  {0,1,-1}
},
{
  {-1,1,2},
  {-1,-1,-1},
  {0,-1,2},
  {0,1,-1}
}
};


template<>
constexpr int8_t t8_standalone_lut_type_face_to_last_facechilds_cubeid<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  {3, 3, 1},
  {3, 3, 2}
};

template<>
constexpr int8_t t8_standalone_lut_type_face_facechildid_to_childid<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACES[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_NUM_FACECHILDREN[T8_ECLASS_TRIANGLE]]={
{
  {1,3},
  {0,3},
  {0,1}
},
{
  {2,3},
  {0,3},
  {0,2}
}
};

template<>
constexpr int8_t t8_standalone_lut_rootface_dim_to_facedim<T8_ECLASS_TRIANGLE>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]={
  {-1,0},
  {0,0},
  {0,-1}
};
template<>
constexpr int8_t t8_standalone_lut_rootface_eq_to_faceeq<T8_ECLASS_TRIANGLE>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]]={
  {-1},  
  {-1},  
  {-1}
};

template<>
constexpr int8_t t8_standalone_lut_type_rootface_to_face<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  {0,1,2},
  {-1,-1,-1}
};

template<>
constexpr int8_t t8_standalone_lut_cornerface<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_CORNERS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_CORNER_FACES[T8_ECLASS_TRIANGLE]]={
  {
  {1,2},
  {0,2},
  {0,1} 
  },
  {
  {1,2},
  {0,2},
  {0,1}  
  }
};

template<>
constexpr int8_t t8_standalone_lut_facecorner<T8_ECLASS_TRIANGLE>[1<<T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_FACE_CORNERS[T8_ECLASS_TRIANGLE]]={
  {
  {1,2},
  {0,2},
  {0,1} 
  },
  {
  {1,2},
  {0,2},
  {0,1} 
  }
};

template<>
constexpr int8_t t8_standalone_lut_transform_coords<T8_ECLASS_TRIANGLE>[1 << T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]={
  {
    {1, 0},
    {0, 1}
  },
  {
    {0, 1},
    {1, 0}
  }
};

template<>
constexpr int8_t t8_standalone_lut_backtransform_coords<T8_ECLASS_TRIANGLE>[1 << T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]={
  {
    {1, 0},
    {0, 1}
  },
  {
    {0, 1},
    {1, 0}
  }
};

template<>
constexpr t8_eclass_t t8_standalone_lut_rootface_to_eclass<T8_ECLASS_TRIANGLE>[T8_ELEMENT_NUM_FACES[T8_ECLASS_TRIANGLE]]={
  T8_ECLASS_LINE, 
  T8_ECLASS_LINE, 
  T8_ECLASS_LINE
};

template <>
constexpr int8_t t8_standalone_lut_bdy_dim_id_idim_to_elem_idim<T8_ECLASS_TRIANGLE>[T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_BOUNDARIES[T8_ECLASS_TRIANGLE]]
                                                       [T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]=
{
  {//bdy_dim 0
    {-1,-1},
    {-1,-1},
    {-1,-1}
  },
  {
    {1,-1}, 
    {0,-1}, //first value could also be 1, since x=y for that face
    {0,-1} 
  }
};

template <>
constexpr int8_t t8_standalone_lut_bdy_dim_id_elem_idim_to_bdy_idim<T8_ECLASS_TRIANGLE>[T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_BOUNDARIES[T8_ECLASS_TRIANGLE]]
                                                       [T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]=
{
  {//bdy_dim 0
    {-1,-1},
    {-1,-1},
    {-1,-1}
  },
  {
    {-1,0}, 
    {0, 0},
    {0,-1} 
  }
};


template <>
constexpr int8_t t8_standalone_lut_bdy_dim_id_elem_idim_to_1_bdy<T8_ECLASS_TRIANGLE>[T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]][T8_ELEMENT_MAX_BOUNDARIES[T8_ECLASS_TRIANGLE]]
                                                       [T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]
{
  {//bdy_dim 0
    {0,0},
    {1,0},
    {1,1}
  },
  {
    {1,0}, 
    {0, 0},
    {0,0} 
  }
};

template <>
constexpr int t8_standalone_lut_line_boundary_id_type_to_equality<T8_ECLASS_TRIANGLE>[T8_ELEMENT_MAX_BOUNDARIES[T8_ECLASS_TRIANGLE]][T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_TRIANGLE]]={};

template <>
constexpr int t8_standalone_lut_num_boundaries<T8_ECLASS_TRIANGLE>[T8_ELEMENT_DIM[T8_ECLASS_TRIANGLE]]=
{3,3};

// clang-format on
#endif /* T8_SELE_LUT_TRIANGLE_HXX */

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

/** \file t8_standalone_elements.hxx
 * Definition of the element class of the standalone scheme.
 */

#ifndef T8_STANDALONE_ELEMENTS_HXX
#define T8_STANDALONE_ELEMENTS_HXX

#include <t8.h>
#include <array>
#include <bitset>

/** Can we delete this? */
#define t8_standalone_element t8_standalone

/** Dimension of the standalone element types */
constexpr uint8_t T8_ELEMENT_DIM[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 3 };

/** Maximum level of the standalone element types
 * \note The maxlevel is lower than 255 so that we can use \ref t8_scheme::element_get_level (uint8_t)
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

/** Number of face children of the standalone element types */
constexpr uint8_t T8_ELEMENT_MAX_NUM_FACECHILDREN[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 4, 4, 4, 4 };

/** Number of faces per corner of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_CORNER_FACES[T8_ECLASS_COUNT] = { 0, 1, 2, 2, 3, 3, 3, 4 };

/** Number of corners per face of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_FACE_CORNERS[T8_ECLASS_COUNT] = { 0, 2, 2, 2, 4, 3, 4, 4 };

/** Number of equations of the standalone element types */
constexpr uint8_t T8_ELEMENT_NUM_EQUATIONS[T8_ECLASS_COUNT] = { 0, 0, 0, 1, 0, 3, 1, 2 };

/* The following lookup tables (LUTs) are only specialized for element classes with type equations
 * (triangle, tet, prism, pyramid), see t8_standalone_lut/. Unless stated otherwise, the first index
 * is the element type interpreted as integer (type.to_ulong ()). */

/**PARENT CHILD BIJECTION*/
/** Type of a child. Indexed by [parent type][child id (Iloc)]. */
template <t8_eclass TEclass>
constexpr int8_t t8_element_type_Iloc_to_childtype[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                  [T8_ELEMENT_NUM_CHILDREN[TEclass]];
/** Cube id of a child, i.e. its position inside the parent's cube. Indexed by [parent type][child id (Iloc)]. */
template <t8_eclass TEclass>
constexpr int8_t t8_element_type_Iloc_to_childcubeid[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                    [T8_ELEMENT_NUM_CHILDREN[TEclass]];
/** Type of the parent. Indexed by [child type][child cube id]. */
template <t8_eclass TEclass>
constexpr int8_t t8_element_type_cubeid_to_parenttype[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                     [1 << T8_ELEMENT_DIM[TEclass]];
/** Child id (Iloc) of an element relative to its parent. Indexed by [child type][child cube id]. */
template <t8_eclass TEclass>
constexpr int8_t t8_element_type_cubeid_to_Iloc[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][1 << T8_ELEMENT_DIM[TEclass]];

/**TYPE EQUATIONS*/
/** The two coordinate dimensions compared by each type equation. Indexed by [equation][0|1].
 * Bit i of the element type encodes the ordering of the coordinates in these two dimensions. */
template <t8_eclass TEclass>
constexpr int8_t t8_type_edge_equations[T8_ELEMENT_NUM_EQUATIONS[TEclass]][2];

/**VERTEX*/
/** Offset (0 or 1, in units of the element length) of a vertex relative to the anchor node.
 * Indexed by [type][vertex][dim]. */
template <t8_eclass TEclass>
constexpr int8_t t8_type_vertex_dim_to_binary[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_NUM_CORNERS[TEclass]]
                                             [T8_ELEMENT_DIM[TEclass]];

/**FACE*/
/** 1 if the face lies inside the element's cube (neighbor differs only in its type), 0 if it lies on the cube boundary.
 * Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_face_internal[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_NUM_FACES[TEclass]];
/** Type bit that is flipped to obtain the neighbor across an internal face, -1 for cube boundary faces.
 * Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_typebit[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                       [T8_ELEMENT_NUM_FACES[TEclass]];
/** Orientation (+1/-1) of the outward face normal along \ref t8_standalone_lut_type_face_to_facenormal_dim,
 * 0 for internal faces. Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_sign[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                    [T8_ELEMENT_NUM_FACES[TEclass]];
/** Coordinate dimension of the face normal for cube boundary faces, -1 for internal faces. Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_facenormal_dim[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                              [T8_ELEMENT_NUM_FACES[TEclass]];
/** Face number of the face neighbor's face that coincides with the given face. Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_neighface[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                         [T8_ELEMENT_NUM_FACES[TEclass]];
/** Face of the parent on which the given face of the child lies, -1 if it lies inside the parent.
 * Indexed by [child type][child cube id][child face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_cubeid_face_to_parentface[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                                 [1 << T8_ELEMENT_DIM[TEclass]]
                                                                 [T8_ELEMENT_NUM_FACES[TEclass]];
/** 1 if the face lies on the upper boundary (coordinate = anchor + length) of the face normal dimension,
 * 0 if it lies on the lower boundary. Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_is_1_boundary[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                             [T8_ELEMENT_NUM_FACES[TEclass]];
/** Cube id of the last descendant touching the face, used to compute the last face descendant.
 * Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_last_facechilds_cubeid[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                                      [T8_ELEMENT_NUM_FACES[TEclass]];
/** Child id of the i-th child touching the face. Indexed by [parent type][parent face][face child id]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_facechildid_to_childid[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                                   [T8_ELEMENT_NUM_FACES[TEclass]]
                                                                   [T8_ELEMENT_MAX_NUM_FACECHILDREN[TEclass]];
/** Face of the child that lies on the given parent face, -1 if none. Indexed by [parent type][child id][parent face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_childid_face_to_childface[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                                 [T8_ELEMENT_NUM_CHILDREN[TEclass]]
                                                                 [T8_ELEMENT_MAX_NUM_FACECHILDREN[TEclass]];
/** Root tree face on which the element face can lie, -1 if an element of this type
 * cannot touch the tree boundary with this face. Indexed by [type][face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_face_to_tree_face[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                         [T8_ELEMENT_NUM_FACES[TEclass]];
/** Coordinate dimension of the boundary face element corresponding to a dimension of the volume element,
 * -1 if the dimension is dropped. Indexed by [root face][volume dim]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_rootface_dim_to_facedim[T8_ELEMENT_NUM_FACES[TEclass]][T8_ELEMENT_DIM[TEclass]];
/** Type equation of the boundary face element corresponding to a type equation of the volume element,
 * -1 if the equation is dropped. Indexed by [root face][volume equation]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_rootface_eq_to_faceeq[T8_ELEMENT_NUM_FACES[TEclass]]
                                                        [T8_ELEMENT_NUM_EQUATIONS[TEclass]];
/** Element face that lies on the given root tree face, -1 if an element of this type cannot touch it.
 * Indexed by [type][root face]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_type_rootface_to_face[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]]
                                                        [T8_ELEMENT_NUM_FACES[TEclass]];
/** Element class of each root tree face. Indexed by [root face]. */
template <t8_eclass TEclass>
constexpr t8_eclass t8_standalone_lut_rootface_to_eclass[T8_ELEMENT_NUM_FACES[TEclass]];

/** Faces adjacent to a corner. Indexed by [type][corner][i]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_cornerface[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_NUM_CORNERS[TEclass]]
                                             [T8_ELEMENT_NUM_CORNER_FACES[TEclass]];

/** Corners of a face. Indexed by [type][face][i]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_facecorner[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_NUM_FACES[TEclass]]
                                             [T8_ELEMENT_NUM_FACE_CORNERS[TEclass]];
/** Matrix mapping reference coordinates of the reference element to the element of the given type.
 * Indexed by [type][out dim][in dim]. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_transform_coords[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_DIM[TEclass]]
                                                   [T8_ELEMENT_DIM[TEclass]];

/** Matrix mapping the internal coordinate frame back to the reference coordinates of the tree.
 * Indexed by [type][out dim][in dim]; currently only type 0 is used. */
template <t8_eclass TEclass>
constexpr int8_t t8_standalone_lut_backtransform_coords[1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]][T8_ELEMENT_DIM[TEclass]]
                                                       [T8_ELEMENT_DIM[TEclass]];

#include "t8_standalone_lut/t8_standalone_lut_triangle.hxx"
#include "t8_standalone_lut/t8_standalone_lut_prism.hxx"
#include "t8_standalone_lut/t8_standalone_lut_pyra.hxx"
#include "t8_standalone_lut/t8_standalone_lut_tet.hxx"

/** Type for the element coordinates. */
typedef int32_t t8_element_coord;
/** Type for the element level. */
typedef uint8_t t8_element_level;
/** Type for the cube id. */
typedef uint8_t t8_cube_id;
/** Type for the child id. */
typedef uint8_t t8_child_id;

/** Define type for the element type. 
 * \tparam TEclass The shape of the element as an eclass.
 */
template <t8_eclass_t TEclass>
using t8_element_type = std::bitset<T8_ELEMENT_NUM_EQUATIONS[TEclass]>;

/** Define type for the element coordinates. 
 * \tparam TEclass The shape of the element as an eclass.
 */
template <t8_eclass_t TEclass>
using t8_element_coords = std::array<t8_element_coord, T8_ELEMENT_DIM[TEclass]>;

/** The data container describing a refined element in a refined tree, where the root element has class \a TEclass */
template <t8_eclass_t TEclass>
struct t8_standalone_element
{
  /** The coordinates of the anchor vertex of the element. */
  t8_element_coords<TEclass> coords;
  /** The refinement level of the element relative to the root at level 0. */
  t8_element_level level;
  /** Bit array: which inequality is fulfilled at which level. */
  t8_element_type<TEclass> type;
};

#endif /* T8_STANDALONE_ELEMENTS_HXX */

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

/** \file t8_dpyramid_bits.h
 * Definitions of pyramid-specific functions.
 */

#ifndef T8_DPYRAMID_BITS_H
#define T8_DPYRAMID_BITS_H

#include "t8_element.h"
#include "t8_dpyramid.h"

T8_EXTERN_C_BEGIN ();

/** Initialize a pyramid as the pyramid with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] p  Existing pyramid whose data will be filled.
 * \param [in] level  level of uniform grid to be considered.
 * \param [in] id     Index to be considered.
 */
void
t8_dpyramid_init_linear_id (t8_dpyramid_t *p, const int level, t8_linearidx_t id);

/** Compute the level of a pyramid.
 * \param [in] p    Pyramid whose level is computed.
 * \return          The level of \a p.
 */
int
t8_dpyramid_get_level (const t8_dpyramid_t *p);

/** Copy the data from source to dest
 * \param[in] source    The source-pyramid
 * \param[in,out] dest  The destination
 */
void
t8_dpyramid_copy (const t8_dpyramid_t *source, t8_dpyramid_t *dest);

/** Computes the linear position of a pyramid in an uniform grid.
 * \param [in] p          pyramid whose id will be computed.
 * \param [in] level      The level on which the linear-id should be computed.
 * \return                Returns the linear position of this pyramid on a grid.
 */
t8_linearidx_t
t8_dpyramid_linear_id (const t8_dpyramid_t *p, const int level);

/** Compute the child_id-th child in Morton order of a pyramid.
 * \param [in] elem         Input pyramid.
 * \param [in,out] child_id The id of the child, 0..7 in Morton order.
 * \param [out] child       Existing pyramid whose data will be filled
 * 		                      with the date of t's child_id-th child.
 */
void
t8_dpyramid_child (const t8_dpyramid_t *elem, const int child_id, t8_dpyramid_t *child);

/** Compute the children of a pyramid, array version
 * \param [in] p        Input pyramid
 * \param [in, out] c   Pointers to the computed children in Morton order
 */
void
t8_dpyramid_children (const t8_dpyramid_t *p, t8_dpyramid_t **c);

/** Given a pyramid and a face, compute all children touching this face
 * \param [in] p                  Input pyramid
 * \param [in] face               The face to compute the children at
 * \param [in, out] children      The children of \a p at \a face
 * \param [in] num_children       The number of children at this face
 * \param [in, out] child_indices An array to be filled with the local-ids of the children.
 */
void
t8_dpyramid_children_at_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *children[], const int num_children,
                              int *child_indices);

/** Given a face of a pyramid and a child number of a child of that face,
 * return the face number of the child of the pyramid that matches the child
 * face.
 * \param[in] p               Input pyramid
 * \param[in] face            A face of \a p
 * \param [in] face_child     A number specifying a child on the \a face
 * \returns                   The number of the face of the child \a face_child*/
int
t8_dpyramid_face_child_face (const t8_dpyramid_t *p, const int face, const int face_child);

/** Given the facenumber of a pyramid, return the shape of the face
 * \param[in] pyra       Input pyramid
 * \param[in] face       The facenumber
 * \return               the shape of the face  
*/
t8_element_shape_t
t8_dpyramid_face_shape (const t8_dpyramid_t *pyra, int face);

/** Returns the corner number of a pyramid given a face of a pyramid and
 *  a corner number regarding that face.
 * \param [in]  pyra    Input pyramid
 * \param [in]  face    The facenumber of a face of \a pyra
 * \param [in]  corner  The cornernumber of a corner of \a face
 * \return              The cornernumber of \a pyra */
int
t8_dpyramid_get_face_corner (const t8_dpyramid_t *pyra, int face, int corner);

/** Given a boundary element and a facenumber of this element, compute the boundary face
 * \param[in] p                   Input pyramid
 * \param[in] face                The face number of an element
 * \param[in, out] boundary       The boundary face*/
void
t8_dpyramid_boundary_face (const t8_dpyramid_t *p, const int face, t8_element_t *boundary);

/** Given a boundary face inside the root pyramids's face construct the element inside the root pyramid that has the 
 * given face as a face.
 * \param [in] face      A face element.
 * \param [in,out] p     An allocated element. The entries will be filled with the data of the element that has \a face
 *                       as a face and lies within the root tree.
 * \param [in] root_face The index of the face of the root tree in which \a face lies.
 * \return               The face number of the face of \a p that coincides with \a face.
 */
int
t8_dpyramid_extrude_face (const t8_element_t *face, t8_dpyramid_t *p, const int root_face);

/** Compare two elements. returns negative if p1 < p2, zero if p1 equals p2 and positive if p1 > p2. 
 * If p2 is a copy of p1 then the elements are equal.
 * \param[in] p1    A pyramid
 * \param[in] p2    Another pyramid
 * \returns         an integer describing which pyramid is larger.
 */
int
t8_dpyramid_compare (const t8_dpyramid_t *p1, const t8_dpyramid_t *p2);

/** Check if two elements are equal.
* \param [in] elem1  The first element.
* \param [in] elem2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dpyramid_equal (const t8_dpyramid_t *elem1, const t8_dpyramid_t *elem2);

/** Check whether a collection of 10 pyramids is a family in Morton order.
 * \param [in]  fam A collection of pyramids
 * \return      Nonzero if \a fam is a family of pyramids
 */
int
t8_dpyramid_is_family (t8_dpyramid_t **fam);

/** Compute whether a given pyramid shares a given face with its root tree.
 * \param [in]  p       The input pyramid
 * \param [in]  face    A face of \a p
 * \return              True, if \a is a subface of the pyramid root element.
 */
int
t8_dpyramid_is_root_boundary (const t8_dpyramid_t *p, const int face);

/** Compute the neighbor of p along a given face and the number of the dual face if
 * the neighbor is inside the root pyramid. Return 0 if the neighbor is not inside, 1 ow.
 * \param[in] p                 Input pyramid
 * \param[in, out] neigh        The neighbor of \a p
 * \param[in] face              The face of \a p along which \a neigh is computed
 * \param [in, out] neigh_face  The dual face
 * */
int
t8_dpyramid_face_neighbor_inside (const t8_dpyramid_t *p, t8_dpyramid_t *neigh, const int face, int *neigh_face);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] p  pyramid to be considered.
 * \return Returns its child id in 0..9
 */
int
t8_dpyramid_child_id (const t8_dpyramid_t *p);

/** Returns zero if p is not inside root, 1 ow
 \param[in] p       Pyramid to check
 \returns 0 if p is inside root, 1, ow*/
int
t8_dpyramid_is_inside_root (const t8_dpyramid_t *p);

/** Check, if a tet of type 0 or 3 has a common face with its pyramid-ancestor
 * \param [in] p      input pyramid
 * \param [in] face   A face of \a p.
 * \return            false if they don't share a face, true otherwise
 */
int
t8_dpyramid_tet_boundary (const t8_dpyramid_t *p, const int face);

/** compute if a given element lies on the tree boundary and return the face number of the tree face. 
 * If not the return value is arbitrary
 * \param [in] p     pyramid
 * \param [in] face     a face of \a p
 * \return              See description
 */
int
t8_dpyramid_tree_face (const t8_dpyramid_t *p, const int face);

/** Compute the first descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the smallest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] desc    Existing pyramid whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void
t8_dpyramid_first_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc, const int level);

/** Construct the first descendant of a pyramid touching a given face
 * \param [in] p           pyramid whose descendant is computed.
 * \param [in] face        The face at which the descendant is computed
 * \param [out] first_desc Existing pyramid whose data will be filled with the data
 *                         of \a p's first descendant on level \a level.
 * \param [in] level       The refinement level. Must be greater than \a p's refinement
 *                         level.*/
void
t8_dpyramid_first_descendant_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *first_desc, const int level);

/** Compute the last descendant of a pyramid at a given level. This is the descendant of
 * the pyramid in a uniform level refinement that has the largest id.
 * \param [in] p        pyramid whose descendant is computed.
 * \param [out] desc    Existing pyramid whose data will be filled with the data of \a p's last descendant on level 
 *                      \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement level.
 */
void
t8_dpyramid_last_descendant (const t8_dpyramid_t *p, t8_dpyramid_t *desc, int level);

/** Construct the last descendant of a pyramid touching a given face
 * \param [in] p            pyramid whose descendant is computed.
 * \param [in] face         The face at which the descendant is computed
 * \param [out] last_desc   Existing pyramid whose data will be filled with the data
 *                          of \a p's first descendant on level \a level.
 * \param [in] level        The refinement level. Must be greater than \a p's refinement level.*/
void
t8_dpyramid_last_descendant_face (const t8_dpyramid_t *p, const int face, t8_dpyramid_t *last_desc, const int level);

/** Compute the coordinates of a vertex of a pyramid.
 * \param [in] elem    Input pyramid.
 * \param [in] vertex  The number of the vertex.
 * \param [out] coords An array of 3 t8_dpyramid_coord_t that will be filled with the coordinates of the vertex.
 */
void
t8_dpyramid_compute_integer_coords (const t8_dpyramid_t *elem, const int vertex, int coords[]);

/** Compute the parent of a given pyramid
 * \param [in] p        Input pyramid.
 * \param [out] parent  The parent of \a p.
 */
void
t8_dpyramid_parent (const t8_dpyramid_t *p, t8_dpyramid_t *parent);

/**
 * Compute the number of corners of a pyramid. If pyramid has type less than 6, it is actually a tetrahedron.
 * \param [in] p    Input pyramid.
 * \return          The number of corners of p.
 */
int
t8_dpyramid_num_corners (const t8_dpyramid_t *p);

/** Compute the number of children of p
 * \param [in] p    Input pyramid.
 * \return          The number of children of p.
 */
int
t8_dpyramid_num_children (const t8_dpyramid_t *p);

/** Compute the number of siblings of p
 * \param [in] p    Input pyramid
 * \return          The number of siblings of p.
 */
int
t8_dpyramid_num_siblings (const t8_dpyramid_t *p);

/** Return the number of faces of p
 * \param [in] p    Input pyramid
 * \return          The number of faces of p
 */
int
t8_dpyramid_num_faces (const t8_dpyramid_t *p);

/** Return the maximal number of faces of an element p
 * \param [in] p    Input pyramid
 * \return          The maximal number of faces of p*/
int
t8_dpyramid_max_num_faces (const t8_dpyramid_t *p);

/** Given a face of an element return the face number of the parent of the element that matches the element's face. 
 * Or return -1 if no face of the parent matches the face.
 * \param [in] elem Input pyramid
 * \param [in] face a face of \a elem
 * \return          the facenumber of the parent of \a elem matching \a face or -1*/
int
t8_dpyramid_face_parent_face (const t8_dpyramid_t *elem, const int face);

/** Return the child-id of the ancestor of p at level level
 * \param [in] p     Input pyramid
 * \param [in] level The ancestor-level
 * \return           The child-id of the ancestor*/
int
t8_dpyramid_ancestor_id (const t8_dpyramid_t *p, const int level);

/**
 * Compute the ancestor of \a pyra at a given level
 * \param[in] pyra           Input pyramid
 * \param[in] level          Level of the ancestor to compute
 * \param[in, out] ancestor  Allocated element that will be filled with the data of the ancestor.
 */
void
t8_dpyramid_ancestor (const t8_dpyramid_t *pyra, const int level, t8_dpyramid_t *ancestor);

/** Compute the type of a pyramid at a given level. Starting from its own level, we iterate over the levels and 
 * compute the type of this level. If p is a tetrahedron, we compute it in a tetrahedral fashion up unto the last 
 * level where p is a tet and continue in a pyramidal fashion 
 * \param [in] p      Input pyramid
 * \param [in] level  The level at which the type is computed
 * \return            The type of \a p at level \a level. */
int
t8_dpyramid_type_at_level (const t8_dpyramid_t *p, const int level);

/** Returns the shape of the pyramid (pyramid or tetrahedron)
 * \param [in] p    Input pyramid.
 * \return          The eclass of the element
 */
t8_element_shape_t
t8_dpyramid_shape (const t8_dpyramid_t *p);

/** Computes the successor of a pyramid in a uniform grid of level \a level.
 * \param [in] elem  pyramid whose id will be computed.
 * \param [in,out] s Existing pyramid whose data will be filled with the data of \a l's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dpyramid_successor (const t8_dpyramid_t *elem, t8_dpyramid_t *s, const int level);

/** Compute the reference coordinates of a vertex of a pyramid when the tree (level 0 triangle) is embedded in \f$ [0,1]^3 \f$.
 * \param [in] elem    Input pyramid.
 * \param [in] vertex  The number of the vertex.
 * \param [out] coords An array of 3 double that will be filled with the reference coordinates of the vertex.
 */
void
t8_dpyramid_vertex_reference_coords (const t8_dpyramid_t *elem, const int vertex, double coords[]);

/** Convert points in the reference space of a pyramid element to points in the
 *  reference space of the tree (level 0) embedded in \f$ [0,1]^3 \f$.
 * \param [in]  elem       Input pyramid.
 * \param [in]  ref_coords The reference coordinates in the pyramid
 *                         (\a num_coords times \f$ [0,1]^3 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [out] out_coords An array of \a num_coords x 3 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the pyramid.
 */
void
t8_dpyramid_compute_reference_coords (const t8_dpyramid_t *elem, const double *ref_coords, const size_t num_coords,
                                      double *out_coords);

/**
 * Compute the nearest common ancestor of two elements
 * \param [in]      pyra1       The first pyramid
 * \param [in]      pyra2       The second pyramid
 * \param [in,out]  nca         Existing pyramid whose data will be filled with
 *                              the data of \a pyra1 and \a pyra2 nearest common ancestor.
 */
void
t8_dpyramid_nearest_common_ancestor (const t8_dpyramid_t *pyra1, const t8_dpyramid_t *pyra2, t8_dpyramid_t *nca);

/** Query whether all entries of a pyramid are in valid ranges.
 * A pyramid is valid if and only if its triangle and line member are valid.
 * \param [in] p  pyramid to be considered.
 * \return        True, if \a p is a valid pyramid and it is safe to call any function in this file on \a p.
 *                False otherwise.
 */
int
t8_dpyramid_is_valid (const t8_dpyramid_t *p);

T8_EXTERN_C_END ();

#endif /* T8_DPYRAMID_BITS_H */

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

/** \file t8_dprism_bits.h
 * Definition if prism-specific functions.
 */

#ifndef T8_DPRISM_BITS_H
#define T8_DPRISM_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_prism/t8_dprism.h>

T8_EXTERN_C_BEGIN ();

/** Look-up table for the corners of each face. */
extern int t8_dprism_face_corners[5][4];

/** Compute the level of a prism.
 * \param [in] p    Line whose prism is computed.
 * \return          The level of \a p.
 */
int
t8_dprism_get_level (const t8_dprism_t *p);

/** Copy all values from one prism to another.
 * \param [in] p    The prism to be copied.
 * \param [in,out] dest Existing prism whose data will be filled with the data
 *                   of \a p.
 */
void
t8_dprism_copy (const t8_dprism_t *p, t8_dprism_t *dest);

/** Compare two elements. returns negative if p1 < p2, zero if p1 equals p2
 *  and positive if p1 > p2.
 *  If p2 is a copy of p1 then the elements are equal.
 */
int
t8_dprism_compare (const t8_dprism_t *p1, const t8_dprism_t *p2);

/** Check if two elements are equal.
* \param [in] elem1  The first element.
* \param [in] elem2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dprism_equal (const t8_dprism_t *elem1, const t8_dprism_t *elem2);

/** Initialize a prism as the prism with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] p  Existing prism whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void
t8_dprism_init_linear_id (t8_dprism_t *p, int level, t8_linearidx_t id);

/** Computes the successor of a prism in a uniform grid of level \a level.
 * \param [in] p  prism whose id will be computed.
 * \param [in,out] succ Existing prism whose data will be filled with the
 *                data of \a l's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dprism_successor (const t8_dprism_t *p, t8_dprism_t *succ, int level);

/** Compute the parent of a prism.
 * \param [in]  p Input prism.
 * \param [in,out] parent Existing prism whose data will
 *                  be filled with the data of p's parent.
 * \note \a p may point to the same prism as \a parent.
 */
void
t8_dprism_parent (const t8_dprism_t *p, t8_dprism_t *parent);

/** Compute the first descendant of a prism at a given level. This is the descendant of
 * the prism in a uniform level refinement that has the smallest id.
 * \param [in] p        Prism whose descendant is computed.
 * \param [out] desc       Existing prism whose data will be filled with the data
 *                      of \a p's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a p's refinement
 *                      level.
 */
void
t8_dprism_first_descendant (const t8_dprism_t *p, t8_dprism_t *desc, int level);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] p  prism to be considered.
 * \return Returns its child id in 0 - 7
 */
int
t8_dprism_child_id (const t8_dprism_t *p);

/** Check whether a collection of eight prism is a family in Morton order.
 * \param [in] fam  An array of eight prism.
 * \return          Nonzero if \a fam is a family of prism.
 */
int
t8_dprism_is_familypv (t8_dprism_t **fam);

/** Computes the nearest common ancestor of two prisms in the same tree.
 * \param [in]     p1 First input tetrahedron.
 * \param [in]     p2 Second input tetrahedron.
 * \param [in,out] r Existing tetrahedron whose data will be filled.
 * \note \a t1, \a t2, \a r may point to the same tetrahedron.
 */
void
t8_dprism_nearest_common_ancestor (const t8_dprism_t *p1, const t8_dprism_t *p2, t8_dprism_t *r);

/** Constructs the boundary element of a prism at a given face
  * \param [in] p       The input prism.
  * \param [in] face    A face of \a p
  * \param [in, out] boundary The boundary element at \a face of \a p*/
void
t8_dprism_boundary_face (const t8_dprism_t *p, int face, t8_element_t *boundary);

/** Compute whether a given prism shares a given face with its root tree.
 * \param [in] p        The input prism.
 * \param [in] face     A face of \a p.
 * \return              True if \a face is a subface of the prisms's root element.
 */
int
t8_dprism_is_root_boundary (const t8_dprism_t *p, int face);

/** Test if a prism lies inside of the root prism, that is the prism of level 0, anchor node (0,0) and type 0.
 *  \param [in]     p Input prism.
 *  \return true    If \a p lies inside of the root prism.
 */
int
t8_dprism_is_inside_root (t8_dprism_t *p);

/** Compute the childid-th child in Morton order of a prism.
 * \param [in] p          Input prism.
 * \param [in] childid    The id of the child, in 0 - 7, in Morton order.
 * \param [in,out] child  Existing prism whose data will be filled with the date of p's childid-th child.
 */
void
t8_dprism_child (const t8_dprism_t *p, int childid, t8_dprism_t *child);

/** Return the shape of a face.
 * \param [in] face A face id for \a p.
 * \return  The shape of \a face.
 */
t8_element_shape_t
t8_dprism_face_shape (int face);

/** Compute the number of children at a given face.
  * \param [in] face The face number
  * \return     Number of Children at \a face
  */
int
t8_dprism_num_face_children (int face);

/** Compute the face neighbor of a prism.
 * \param [in]     p      Input prism.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] neigh      Existing prism whose data will be filled.
 * \note \a p may point to the same prism as \a n.
 * \return  The face number of \a neigh of the respective face.
 */
int
t8_dprism_face_neighbour (const t8_dprism_t *p, int face, t8_dprism_t *neigh);

/** Return the corner number of a prism corresponding to a given face corner.
  * \param [in] face    The face number.
  * \param [in] corner  A corner of \a face
  * \return             The corner index of \a p corresponding to the \a corner-th corner of \a face.
  */
int
t8_dprism_get_face_corner (int face, int corner);

/** Compute the children of a prism, array version.
 * \param [in] p Pointers to prism element the children should be computed for.
 * \param [in] length Number of children.
 * \param [in,out] c  Pointers to the computed children in Morton order.
 *
 */
void
t8_dprism_childrenpv (const t8_dprism_t *p, int length, t8_dprism_t *c[]);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] p  prism to be considered.
 * \param [in] level level to be considered.
 * \return Returns its child id in 0..7
 */
int
t8_dprism_ancestor_id (t8_dprism_t *p, int level);

/** Given a prism and a face of the prism, compute all children of the prism that touch the face.
 * \param [in] p            The prism.
 * \param [in] face         A face of \a p.
 * \param [in,out] children Allocated prism, in which the children of \a p that share a face with \a face are stored.
 *                          They will be stored in order of their child_id.
 * \param [in] num_children The number of prisms in \a children. Must match the number of children that touch \a face.
 * \param [in,out] child_indices The indices of the children in \a children. Only filled if this is null previously.
 */
void
t8_dprism_children_at_face (const t8_dprism_t *p, int face, t8_dprism_t **children, int num_children,
                            int *child_indices);

/** Given a face of a prism and a child number of a child of that face, return the face number of the child of the
 * prism that matches the child face.
 * \param [in]  face        The number of the face.
 * \return                  The face number of the face of a child of \a p that coincides with \a face_child.
 */
int
t8_dprism_face_child_face (int face);

/** Given a face of a prism return the face number of the parent of the prism that matches the prism's face.
 * Or return -1 if no face of the parent matches the face.
 * \param [in]  prism   The prism.
 * \param [in]  face    Then number of the face.
 * \return              If \a face of \a prism is also a face of \a prism's parent, the face number of this face.
 *                      Otherwise -1.
 */
int
t8_dprism_face_parent_face (const t8_dprism_t *prism, int face);

/** Given a prism and a face of this prism. If the face lies on the tree boundary, return the face number of the tree
 * face. If not the return value is arbitrary.
 * \param [in] face The index of a face of \a elem.
 * \return          The index of the tree face that \a face is a subface of, if \a face is on a tree boundary.
 *                  Any arbitrary integer if \a is not at a tree boundary.
 */
int
t8_dprism_tree_face (int face);

/** Given a boundary face inside a root tree's face construct the element inside the root tree that has the given face
 * as a face.
 * \param [in] face       A face element.
 * \param [in,out] elem   An allocated element. The entries will be filled with
 *                        the data of the element that has \a face as a face and
 *                        lies within the root tree.
 * \param [in] root_face  The index of the face of the root tree in which \a face lies.
 */
void
t8_dprism_extrude_face (const t8_element_t *face, t8_element_t *elem, const int root_face);

/** Compute the last descendant of a prism at a given level. This is the descendant of the prism in a uniform level
 * refinement that has the largest id.
 * \param [in] p      Prism whose descendant is computed.
 * \param [out] s     Existing prism whose data will be filled with the data
 *                    of \a p's last descendant on level \a level.
 * \param [in] level  The refinement level. Must be greater than \a p's refinement level.
 */
void
t8_dprism_last_descendant (const t8_dprism_t *p, t8_dprism_t *s, int level);

/** Compute the descendant of a prism in a given corner.
 * \param [in] p        Prism whose descendant is computed.
 * \param [out] s       Existing prism whose data will be filled with the data
 *                      of p's descendant in \a corner.
 * \param [in]  corner  The corner in which the descendant should lie.
 * \param [in]  level   The refinement level of the descendant. Must be greater or
 *                      equal to \a p's level.
 */
void
t8_dprism_corner_descendant (const t8_dprism_t *p, t8_dprism_t *s, int corner, int level);

/** Computes the linear position of a prism in an uniform grid.
 * \param [in] p  Prism whose id will be computed.
 * \param [in] level  Refinement level of the prism.
 * \return        Returns the linear position of this prism on a grid.
 */
t8_linearidx_t
t8_dprism_linear_id (const t8_dprism_t *p, int level);

/** Query whether all entries of a prism are in valid ranges. A prism is valid if and only if its triangle and
 * line member are valid.
 * \param [in] p  prism to be considered.
 * \return        True, if \a p is a valid prism and it is safe to call any function in this file on \a p.
 *                False otherwise.
 */
int
t8_dprism_is_valid (const t8_dprism_t *p);

T8_EXTERN_C_END ();

#endif /* T8_DPRISM_BITS_H */

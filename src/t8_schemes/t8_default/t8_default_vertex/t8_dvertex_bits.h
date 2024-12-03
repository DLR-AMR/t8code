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

/** \file t8_dvertex_bits.h
 * Definitions of vertex-specific functions.
 */

#ifndef T8_DVERTEX_BITS_H
#define T8_DVERTEX_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_vertex/t8_dvertex.h>

T8_EXTERN_C_BEGIN ();

/** Compute the level of a vertex.
 * \param [in] l    vertex whose level is computed.
 * \return          The level of \a l.
 */
int
t8_dvertex_get_level (const t8_dvertex_t *v);

/** Copy all values from one vertex to another.
 * \param [in] l    The vertex to be copied.
 * \param [in,out] dest Existing vertex whose data will be filled with the data
 *                   of \a l.
 */
void
t8_dvertex_copy (const t8_dvertex_t *v, t8_dvertex_t *dest);

/** Compare two elements. returns negative if l1 < l2, zero if l1 equals l2
 *  and positive if l1 > l2.
 *  If l2 is a copy of l1 then the elements are equal.
 */
int
t8_dvertex_compare (const t8_dvertex_t *l1, const t8_dvertex_t *l2);

/** Check if two elements are equal.
* \param [in] elem1  The first element.
* \param [in] elem2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dvertex_equal (const t8_dvertex_t *elem1, const t8_dvertex_t *elem2);

/** Compute the parent of a vertex.
 * \param [in]  l   The input vertex.
 * \param [in,out] parent Existing vertex whose data will be filled with the parent
 *                  data of \a l.
 */
void
t8_dvertex_parent (const t8_dvertex_t *v, t8_dvertex_t *parent);

/** Compute the childid-th child in Morton order of a vertex.
 * \param [in] l    Input vertex.
 * \param [in] childid The id of the child, 0 or 1, in Morton order.
 * \param [in,out] child  Existing vertex whose data will be filled
 * 		    with the date of l's childid-th child.
 */
void
t8_dvertex_child (const t8_dvertex_t *v, t8_dvertex_t *child);

/** Computes the nearest common ancestor of two vertexs in the same tree.
 * \param [in]     l1 First input vertex.
 * \param [in]     l2 Second input vertex.
 * \param [in,out] r Existing vertex whose data will be filled.
 * \note \a l1, \a l2, \a r may point to the same vertex.
 */
void
t8_dvertex_nearest_common_ancestor (const t8_dvertex_t *t1, const t8_dvertex_t *t2, t8_dvertex_t *r);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] l  vertex to be considered.
 * \param [in] level level to be considered.
 * \return Returns its child id 0 or 1.
 */
int
t8_dvertex_ancestor_id (const t8_dvertex_t *v, int level);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \param [in] t  vertex to be considered.
 * \return Returns its child id in 0,1
 */
int
t8_dvertex_child_id (const t8_dvertex_t *t);

/** Compute the sibling of a vertex.
 * \param [in] v  vertex to be considered.
 * \param [in] sibid The id of the sibling, must be 0.
 * \param [out] s  The sibling of \a v. For vertices \a s is just a copy of \a v.
 */
void
t8_dvertex_sibling (const t8_dvertex_t *v, int sibid, t8_dvertex_t *s);

/** Compute the 2 children of a vertex, array version.
 * \param [in]     t  Input vertex.
 * \param [in,out] c  Pointers to the 2 computed children in Morton order.
 *                    t may point to the same quadrant as c[0].
 */
void
t8_dvertex_childrenpv (const t8_dvertex_t *t, t8_dvertex_t *c[T8_DVERTEX_CHILDREN]);

/** Check whether a collection of two vertexs is a family in Morton order.
 * \param [in]     f  An array of two vertexs.
 * \return            Nonzero if \a f is a family of vertexs.
 */
int
t8_dvertex_is_familypv (const t8_dvertex_t *f[]);

/** Compute whether a given vertex shares a given face with its root tree.
 * \param [in] p        The input vertex.
 * \param [in] face     A face of \a p.
 * \return              True if \a face is a subface of the vertex's root element.
 */
int
t8_dvertex_is_root_boundary (const t8_dvertex_t *p, int face);

/** Test if a vertex lies inside of the root vertex,
 *  that is the vertex of level 0, anchor node (0,0)
 *  \param [in]     l Input vertex.
 *  \return true    If \a l lies inside of the root vertex.
 */
int
t8_dvertex_is_inside_root (const t8_dvertex_t *v);

/** Initialize a vertex as the vertex with a given global id in a uniform
 *  refinement of a given level. *
 * \param [in,out] l  Existing vertex whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void
t8_dvertex_init_linear_id (t8_dvertex_t *v, int level, t8_linearidx_t id);

/** Suppose we have two trees that share a common vertex face f.
 *  Given a vertex e that is a subface of f in one of the trees
 *  and given the orientation of the tree connection, construct the face
 *  vertex of the respective tree neighbor that logically coincides with e
 *  but lies in the coordinate system of the neighbor tree.
 *  \param [in] elem1     The face element.
   *  \param [in,out] elem2 On return the face element \a elem1 with respect
 *                        to the coordinate system of the other tree.
 * \note For vertices this function is equivalent to copy.
 */
void
t8_dvertex_transform_face (const t8_dvertex_t *vertex1, t8_dvertex_t *vertex2);

/** Compute the first descendant of a vertex at a given level. This is the descendant of
 * the vertex in a uniform level refinement that has the smallest id.
 * \param [in] l        vertex whose descendant is computed.
 * \param [out] s       Existing vertex whose data will be filled with the data
 *                      of \a l's first descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a l's refinement
 *                      level.
 */
void
t8_dvertex_first_descendant (const t8_dvertex_t *v, t8_dvertex_t *s, int level);

/** Compute the last descendant of a vertex at a given level. This is the descendant of
 * the vertex in a uniform level refinement that has the largest id.
 * \param [in] l        vertex whose descendant is computed.
 * \param [out] s       Existing vertex whose data will be filled with the data
 *                      of \a l's last descendant on level \a level.
 * \param [in] level    The refinement level. Must be greater than \a l's refinement
 *                      level.
 */
void
t8_dvertex_last_descendant (const t8_dvertex_t *v, t8_dvertex_t *s, int level);

/** Compute the coordinates of a vertex (always 0).
 * \param [in] elem     vertex whose vertex is computed.
 * \param [in] vertex   The number of the vertex of \a elem
 * \param [out] coords   The coordinates of the computed vertex
 */
void
t8_dvertex_vertex_integer_coords (const t8_dvertex_t *elem, int vertex, int coords[]);

/** Compute the coordinates of a vertex (always 0) inside the [0,1]^0 reference space.
 * \param [in] elem     vertex whose vertex is computed.
 * \param [in] vertex   The number of the vertex of \a elem (must be 0).
 * \param [out] coords  The coordinates of the computed vertex, must have one entry (will be set to 0).
 */
void
t8_dvertex_vertex_ref_coords (const t8_dvertex_t *elem, int vertex, double coords[]);

/** Convert points in the reference space of a vertex element to points in the
 *  reference space of the tree (level 0) embedded in \f$ [0,1]^1 \f$.
 * \param [in]  elem       Input vertex.
 * \param [in]  ref_coords The reference coordinates in the vertex
 *                         (\a num_coords times \f$ [0,1]^1 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [in]  padding    Only used if \a num_coords > 1.
 *                         The number of padding entries in \a ref_coords between array elements:
 *                         For elem dim 2 and input {x1, y1, x2, y2, ...} padding is 0.
 *                         For elem dim 2 and input {x1, y1, z1, x2, y2, z2, ...} padding is 1.
 *                         For elem dim 3 and input {x1, y1, z1, x2, y2, z2, ...} padding is 0.
 * \param [out] out_coords An array of \a num_coords x 1 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the vertex (will be set to 0).
 */
void
t8_dvertex_compute_reference_coords (const t8_dvertex_t *elem, const double *ref_coords, const size_t num_coords,
                                     const size_t padding, double *out_coords);

/** Computes the linear position of a vertex in an uniform grid.
 * \param [in] vertex  vertex whose id will be computed.
 * \return Returns the linear position of this vertex on a grid.
 */
t8_linearidx_t
t8_dvertex_linear_id (const t8_dvertex_t *elem, int level);

/** Query whether all entries of a vertex are in valid ranges.
 * \param [in] l  vertex to be considered.
 * \return        True, if \a l is a valid vertex and it is safe to call any
 *                function in this file on \a l.
 *                False otherwise.
 */
int
t8_dvertex_is_valid (const t8_dvertex_t *v);

/** Print a vertex
 * \param [in] v  vertex to be considered.
 */
void
t8_dvertex_debug_print (const t8_dvertex_t *v);

/** Set default values for a vertex, such that it passes \ref t8_dvertex_is_valid.
 * \param [in] l  vertex to be initialized
 */
void
t8_dvertex_init (t8_dvertex_t *v);

T8_EXTERN_C_END ();

#endif /* T8_DVERTEX_BITS_H */

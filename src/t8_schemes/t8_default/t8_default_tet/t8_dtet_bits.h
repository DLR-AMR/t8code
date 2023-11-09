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

/** \file t8_dtet_bits.h
 * Definitions of tet-specific functions.
 * TODO: Run make doxygen and grep for files.
 *       Also document all arguments of functions.
 * TODO: Group the dtet_is functions together.
 */

#ifndef T8_DTET_BITS_H
#define T8_DTET_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_tet/t8_dtet.h>

T8_EXTERN_C_BEGIN ();

/** Compute the coordinates of a vertex of a tetrahedron.
 * \param [in] t    Input tetrahedron.
 * \param [in] vertex The number of the vertex.
 * \param [out] coordinates An array of 3 t8_dtet_coord_t that will be filled with the coordinates of the vertex.
 */
void
t8_dtet_compute_integer_coords (const t8_dtet_t *elem, int vertex, t8_dtet_coord_t coordinates[3]);

/** Compute the coordinates of a vertex of a tetrahedron when the 
 * tree (level 0 tetrahedron) is embedded in \f$ [0,1]^3 \f$.
 * \param [in] elem         Input tetrahedron.
 * \param [in] vertex       The number of the vertex.
 * \param [out] coordinates An array of 3 double that will be filled with the reference coordinates of the vertex.
 */
void
t8_dtet_compute_vertex_ref_coords (const t8_dtet_t *elem, int vertex, double coordinates[3]);

/** Convert points in the reference space of a tet element to points in the
 *  reference space of the tree (level 0) embedded in \f$ [0,1]^3 \f$.
 * \param [in]  elem       Input tet.
 * \param [in]  ref_coords The reference coordinates in the tet
 *                         (\a num_coords times \f$ [0,1]^3 \f$)
 * \param [in]  num_coords Number of coordinates to evaluate
 * \param [out] out_coords An array of \a num_coords x 3 x double that
 * 		                     will be filled with the reference coordinates
 *                         of the points on the tet.
 */
void
t8_dtet_compute_reference_coords (const t8_dtet_t *elem, const double *ref_coords, const size_t num_coords,
                                  double *out_coords);

/** Compute the coordinates of the four vertices of a tetrahedron.
 * \param [in] elem         Input tetrahedron.
 * \param [out] coordinates An array of 4x3 t8_dtet_coord_t that will be filled with the coordinates of t's vertices.
 */
void
t8_dtet_compute_all_coords (const t8_dtet_t *elem, t8_dtet_coord_t coordinates[4][3]);

/** Copy the values of one tetrahedron to another.
 * \param [in] t Tetrahedron whose values will be copied.
 * \param [in,out] dest Existing tetrahedron whose data will be filled with the data of \a t. *
 */
void
t8_dtet_copy (const t8_dtet_t *t, t8_dtet_t *dest);

/** Compare two tets in their linear order.
 * \param [in] t1 Tetrahedron one.
 * \param [in] t2 Tetrahedron two.
 * \return        Returns negative if t1 < t2, zero if t1 = t2, positive if t1 > t2
 */
int
t8_dtet_compare (const t8_dtet_t *t1, const t8_dtet_t *t2);

/** Check if two elements are equal.
* \param [in] elem1  The first element.
* \param [in] elem2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dtet_equal (const t8_dtet_t *elem1, const t8_dtet_t *elem2);

/** Compute the parent of a tetrahedron.
 * \param [in]  elem Input tetrahedron.
 * \param [in,out] parent Existing tetrahedron whose data will be filled with the data of elem's parent.
 * \note \a elem may point to the same tetrahedron as \a parent.
 */
void
t8_dtet_parent (const t8_dtet_t *t, t8_dtet_t *parent);

/** Compute the ancestor of a tetrahedron at a given level.
 * \param [in]  t   Input tetrahedron.
 * \param [in]  level A smaller level than \a t.
 * \param [in,out] ancestor Existing tetrahedron whose data will be filled with the data of \a t's ancestor on
 *                  level \a level.
 * \note The tetrahedron \a ancestor may point to the same tetrahedron as \a t.
 */
void
t8_dtet_ancestor (const t8_dtet_t *t, int level, t8_dtet_t *ancestor);

/** Compute the childid-th child in Morton order of a tetrahedron t.
 * \param [in] t    Input tetrahedron.
 * \param [in,out] childid The id of the child, 0..7 in Bey order.
 * \param [out] child  Existing tetrahedron whose data will be filled with the date of t's childid-th child.
 */
void
t8_dtet_child (const t8_dtet_t *elem, int childid, t8_dtet_t *child);

/** Compute the 8 children of a tetrahedron, array version.
 * \param [in]     t  Input tetrahedron.
 * \param [in,out] c  Pointers to the 8 computed children in Morton order. t may point to the same quadrant as c[0].
 */
void
t8_dtet_childrenpv (const t8_dtet_t *t, t8_dtet_t *c[T8_DTET_CHILDREN]);

/** Check whether a collection of eight tetrahedra is a family in Morton order.
 * \param [in]     f  An array of eight tetrahedra.
 * \return            Nonzero if \a f is a family of tetrahedra.
 */
int
t8_dtet_is_familypv (const t8_dtet_t *f[]);

/** Compute a specific sibling of a tetrahedron.
 * \param [in]     elem  Input tetrahedron.
 * \param [in,out] sibling  Existing tetrahedron whose data will be filled with the data of sibling no. sibling_id of 
 *                          elem.
 * \param [in]     sibid The id of the sibling computed, 0..7 in Bey order.
 */
void
t8_dtet_sibling (const t8_dtet_t *elem, int sibid, t8_dtet_t *sibling);

/** Compute the face neighbor of a tetrahedron.
 * \param [in]     t      Input tetrahedron.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] n      Existing tetrahedron whose data will be filled.
 * \note \a t may point to the same tetrahedron as \a n.
 */
int
t8_dtet_face_neighbour (const t8_dtet_t *t, int face, t8_dtet_t *n);

/** Computes the nearest common ancestor of two tetrahedra in the same tree.
 * \param [in]     t1 First input tetrahedron.
 * \param [in]     t2 Second input tetrahedron.
 * \param [in,out] r Existing tetrahedron whose data will be filled.
 * \note \a t1, \a t2, \a r may point to the same tetrahedron.
 */
void
t8_dtet_nearest_common_ancestor (const t8_dtet_t *t1, const t8_dtet_t *t2, t8_dtet_t *r);

/** Given a tetrahedron and a face of the tetrahedron, compute all children of the tetrahedron that touch the face.
 * \param [in] tet      The tetrahedron.
 * \param [in] face     A face of \a tet.
 * \param [in,out] children Allocated tetrahedra, in which the children of \a tet that share a face with \a face are
 *                      stored. They will be stored in order of their child_id.
 * \param [in] num_children The number of tetrahedra in \a children. Must match the number of children that touch 
 *                      \a face.
 */
void
t8_dtet_children_at_face (const t8_dtet_t *tet, int face, t8_dtet_t *children[], int num_children, int *child_indices);

/** Given a face of an tetrahedron and a child number of a child of that face, return the face number
 * of the child of the tetrahedron that matches the child face.
 * \param [in]  tet     The tetrahedron.
 * \param [in]  face    Then number of the face.
 * \param [in]  face_child  The child number of a child of the face tetrahedron.
 * \return              The face number of the face of a child of \a tetrahedron that coincides with \a face_child.
 */
int
t8_dtet_face_child_face (const t8_dtet_t *tet, int face, int face_child);

/** Given a face of an tet return the face number of the parent of the tet that matches the tet's face. Or return -1 if
 * no face of the parent matches the face.
 * \param [in]  elem  The tet.
 * \param [in]  face  Then number of the face.
 * \return            If \a face of \a elem is also a face of \a elem's parent, the face number of this face. 
 *                    Otherwise -1.
 */
int
t8_dtet_face_parent_face (const t8_dtet_t *tet, int face);

/** Given a tetrahedron and a face of this tetrahedron. If the face lies on the tree boundary, return the face number 
 * of the tree face. If not the return value is arbitrary.
 * \param [in] t    The tetrahedron.
 * \param [in] face The index of a face of \a t.
 * \return The index of the tree face that \a face is a subface of, if \a face is on a tree boundary.
 *         Any arbitrary integer if \a is not at a tree boundary.
 * \note For boundary tetrahedra, this function is the inverse of \ref t8_dtet_root_face_to_face.
 */
int
t8_dtet_tree_face (t8_dtet_t *t, int face);

/** Given a tetrahedron and a face of the root tetrahedron. If the tetrahedron lies on the tree boundary, 
 * return the corresponding face number of the tetrahedron. If not the return value is arbitrary.
 * \param [in] t    The tetrahedron.
 * \param [in] face The index of a face of the root tetrahedron.
 * \return The index of the face of \a t that is a subface of \a face, if \a t is on the tree boundary.
 *         Any arbitrary integer if \a t is not at a tree boundary.
 * \note For boundary tetrahedra, this function is the inverse of \ref t8_dtet_tree_face.
 */
int
t8_dtet_root_face_to_face (t8_dtet_t *t, int root_face);

/** Test if a tetrahedron lies inside of the root tetrahedron, that is the tetrahedron of level 0, anchor node (0,0,0)
 *  and type 0.
 *  \param [in] t   Input tetrahedron.
 *  \return true    If \a t lies inside of the root tetrahedron.
 */
int
t8_dtet_is_inside_root (t8_dtet_t *t);

/** Compute whether a given tetrahedron shares a given face with its root tree.
 * \param [in] t    The input tet.
 * \param [in] face A face of \a t.
 * \return          True if \a face is a subface of the tet's root element.
 */
int
t8_dtet_is_root_boundary (const t8_dtet_t *t, int face);

/** Test if two tetrahedra have the same coordinates, type and level.
 * \return true if \a t1 describes the same tetrahedron as \a t2.
 */
int
t8_dtet_is_equal (const t8_dtet_t *t1, const t8_dtet_t *t2);

/** Test if two tetrahedra are siblings.
 * \param [in] t1 First tetrahedron to be tested.
 * \param [in] t2 Second tetrahedron to be tested.
 * \return true if \a t1 is equal to or a sibling of \a t2.
 */
int
t8_dtet_is_sibling (const t8_dtet_t *t1, const t8_dtet_t *t2);

/** Test if a tetrahedron is the parent of another tetrahedron.
 * \param [in] t tetrahedron to be tested.
 * \param [in] c Possible child tetrahedron.
 * \return true if \a t is the parent of \a c.
 */
int
t8_dtet_is_parent (const t8_dtet_t *t, const t8_dtet_t *c);

/** Test if a tetrahedron is an ancestor of another tetrahedron.
 * \param [in] t tetrahedron to be tested.
 * \param [in] c Descendent tetrahedron.
 * \return true if \a t is equal to or an ancestor of \a c.
 */
int
t8_dtet_is_ancestor (const t8_dtet_t *t, const t8_dtet_t *c);

/** Computes the linear position of a tetrahedron in a uniform grid.
 * \param [in] t  tetrahedron whose id will be computed.
 * \param [in] level level of uniform grid to be considered.
 * \return Returns the linear position of this tetrahedron on a grid of level \a level.
 * \note This id is not the Morton index.
 */
t8_linearidx_t
t8_dtet_linear_id (const t8_dtet_t *t, int level);

/**
 * Same as init_linear_id, but we only consider the subtree. Used for computing the index of a
 * tetrahedron lying in a pyramid
 * \param [in, out] t   Existing triangle whose data will be filled
 * \param id            Index to be considered
 * \param start_level   The level of the root of the subtree
 * \param end_level     Level of uniform grid to be considered
 * \param parenttype    The type of the parent.
 */
void
t8_dtet_init_linear_id_with_level (t8_dtet_t *t, t8_linearidx_t id, int start_level, int end_level,
                                   t8_dtet_type_t parenttype);

/** Initialize a tetrahedron as the tetrahedron with a given global id in a uniform refinement of a given level.
 * \param [in,out] t  Existing tetrahedron whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void
t8_dtet_init_linear_id (t8_dtet_t *t, t8_linearidx_t id, int level);

/** Initialize a tetrahedron as the root tetrahedron (type 0 at level 0)
 * \param [in,out] t Existing tetrahedron whose data will be filled.
 */
void
t8_dtet_init_root (t8_dtet_t *t);

/** Computes the successor of a tetrahedron in a uniform grid of level \a level.
 * \param [in] t  tetrahedron whose id will be computed.
 * \param [out] s Existing tetrahedron whose data will be filled with the
 *                data of t's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dtet_successor (const t8_dtet_t *t, t8_dtet_t *s, int level);

/** Compute the first descendant of a tetrahedron at a given level. This is the descendant of
 * the tetrahedron in a uniform maxlevel refinement that has the smaller id.
 * \param [in] t        tetrahedron whose descendant is computed.
 * \param [in] level    A given level. Must be greater or equal to \a t's level.
 * \param [out] s       Existing tetrahedron whose data will be filled with the data
 *                      of t's first descendant.
 */
void
t8_dtet_first_descendant (const t8_dtet_t *t, t8_dtet_t *s, int level);

/** Compute the last descendant of a tetrahedron at a given level. This is the descendant of
 * the tetrahedron in a uniform maxlevel refinement that has the bigges id.
 * \param [in] t        tetrahedron whose descendant is computed.
 * \param [in] level    A given level. Must be greater or equal to \a t's level.
 * \param [out] s       Existing tetrahedron whose data will be filled with the data of t's last descendant.
 */
void
t8_dtet_last_descendant (const t8_dtet_t *t, t8_dtet_t *s, int level);

/** Compute the descendant of a tetrahedron in a given corner.
 * \param [in] t        Tetrahedron whose descendant is computed.
 * \param [out] s       Existing tetrahedron whose data will be filled with the data of t's descendant in \a corner.
 * \param [in]  corner  The corner in which the descendant should lie.
 * \param [in]  level   The refinement level of the descendant. Must be greater or equal to \a t's level.
 */
void
t8_dtet_corner_descendant (const t8_dtet_t *t, t8_dtet_t *s, int corner, int level);

/** Computes the predecessor of a tetrahedron in a uniform grid of level \a level.
 * \param [in] t     tetrahedron whose id will be computed.
 * \param [in,out] s Existing tetrahedron whose data will be filled with the data of t's predecessor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dtet_predecessor (const t8_dtet_t *t, t8_dtet_t *s, int level);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] t  tetrahedron to be considered.
 * \param [in] level level to be considered.
 * \return Returns its child id in 0..7
 */
int
t8_dtet_ancestor_id (const t8_dtet_t *t, int level);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] t  tetrahedron to be considered.
 * \return Returns its child id in 0..7
 */
int
t8_dtet_child_id (const t8_dtet_t *t);

/** Return the level of a tetrahedron.
 * \param [in] t  tetrahedron to be considered.
 * \return        The level of \a t.
 */
int
t8_dtet_get_level (const t8_dtet_t *t);

/** Query whether all entries of a tet are in valid ranges.
 * \param [in] t  tet to be considered.
 * \return        True, if \a t is a valid tet and it is safe to call any function on \a t. False otherwise.
 */
int
t8_dtet_is_valid (const t8_dtet_t *t);

/** Set sensible default values for a tet.
 * \param [in,out] t A tet.
 */
void
t8_dtet_init (t8_dtet_t *t);

void
t8_dtet_element_pack (t8_dtet_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                      int *position, sc_MPI_Comm comm);

void
t8_dtet_element_pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size);

void
t8_dtet_element_unpack (void *recvbuf, const int buffer_size, int *position, t8_dtet_t **elements,
                        const unsigned int count, sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* T8_DTET_BITS_H */

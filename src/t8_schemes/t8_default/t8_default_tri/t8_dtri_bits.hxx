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

/** \file t8_dtri_bits.hxx
 * Definitions of triangle-specific functions.
 */

#ifndef T8_DTRI_BITS_H
#define T8_DTRI_BITS_H

#include <t8_element.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.hxx>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>

T8_EXTERN_C_BEGIN ();

/** Copy the values of one triangle to another.
 * \param [in] element Triangle whose values will be copied.
 * \param [in,out] dest Existing triangle whose data will be
 *                      filled with the data of \a element.
 */
void
t8_dtri_copy (const t8_dtri_t *element, t8_dtri_t *dest);

/** Compare two triangle in their linear order.
 * \param [in] element1 Triangle one.
 * \param [in] element2 Triangle two.
 * \return        Returns negative if tri1 < tri2, zero if tri1 = tri2, positive if tri1 > tri2
 */
int
t8_dtri_compare (const t8_dtri_t *element1, const t8_dtri_t *element2);

/** Check if two elements are equal.
* \param [in] element1  The first element.
* \param [in] element2  The second element.
* \return            1 if the elements are equal, 0 if they are not equal
*/
int
t8_dtri_equal (const t8_dtri_t *element1, const t8_dtri_t *element2);

/** Compute the parent of a triangle.
 * \param [in]  element Input triangle.
 * \param [in,out] parent Existing triangle whose data will be filled with the data of elem's parent.
 * \note \a element may point to the same triangle as \a parent.
 */
void
t8_dtri_parent (const t8_dtri_t *element, t8_dtri_t *parent);

/** Compute the ancestor of a triangle at a given level.
 * \param [in]  element   Input triangle.
 * \param [in]  level A smaller level than \a element.
 * \param [in,out] ancestor Existing triangle whose data will be filled with the data of \a element's ancestor on level
 *                 \a level.
 * \note The triangle \a ancestor may point to the same triangle as \a element.
 */
void
t8_dtri_ancestor (const t8_dtri_t *element, int level, t8_dtri_t *ancestor);

/** Compute the childid-th child in Morton order of a triangle.
 * \param [in] element    Input triangle.
 * \param [in,out] childid The id of the child, 0..7 in Morton order.
 * \param [out] child  Existing triangle whose data will be filled with the date of t's childid-th child.
 */
void
t8_dtri_child (const t8_dtri_t *element, int childid, t8_dtri_t *child);

/** Compute the 4 children of a triangle, array version.
 * \param [in]     element  Input triangle.
 * \param [in,out] children  Pointers to the 4 computed children in Morton order. t may point to the same quadrant as c[0].
 */
void
t8_dtri_childrenpv (const t8_dtri_t *element, t8_dtri_t *children[T8_DTRI_CHILDREN]);

/** Check whether a collection of eight triangles is a family in Morton order.
 * \param [in]     f  An array of eight triangles.
 * \return            Nonzero if \a f is a family of triangles.
 */
int
t8_dtri_is_familypv (const t8_dtri_t *f[]);

/** Compute a specific sibling of a triangle.
 * \param [in]     element  Input triangle.
 * \param [in]     sibid The id of the sibling computed, 0..7 in Bey order.
 * \param [in,out] sibling  Existing triangle whose data will be filled with the data of sibling no. sibling_id of
 *                          \a element.
 */
void
t8_dtri_sibling (const t8_dtri_t *element, int sibid, t8_dtri_t *sibling);

/** Compute the face neighbor of a triangle.
 * \param [in]     element      Input triangle.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] neigh      Existing triangle whose data will be filled.
 * \note \a element may point to the same triangle as \a neigh.
 */
int
t8_dtri_face_neighbour (const t8_dtri_t *element, int face, t8_dtri_t *neigh);

/** Computes the nearest common ancestor of two triangles in the same tree.
 * \param [in]     element1 First input triangle.
 * \param [in]     element2 Second input triangle.
 * \param [in,out] nca Existing triangle whose data will be filled.
 * \note \a element1, \a element2, \a nca may point to the same quadrant.
 */
void
t8_dtri_nearest_common_ancestor (const t8_dtri_t *element1, const t8_dtri_t *element2, t8_dtri_t *nca);

/** Given a triangle and a face of the triangle, compute all children of the triangle that touch the face.
 * \param [in] element      The triangle.
 * \param [in] face     A face of \a element.
 * \param [in,out] children Allocated triangles, in which the children of \a element that share a face with \a face are
 *                      stored. They will be stored in order of their child_id.
 * \param [in] num_children The number of triangles in \a children. Must match the number of children that touch
 *                      \a face.
 * \param [in,out] child_indices The indices of the children in \a children. Only filled if this is null previously.
 */
void
t8_dtri_children_at_face (const t8_dtri_t *element, int face, t8_dtri_t *children[], int num_children,
                          int *child_indices);

/** Given a face of a triangle and a child number of a child of that face, return the face number of the child of the
 * triangle that matches the child face.
 * \param [in]  triangle The triangle.
 * \param [in]  face    Then number of the face.
 * \param [in]  face_child  The child number of a child of the face triangle.
 * \return              The face number of the face of a child of \a triangle that coincides with \a face_child.
 */
int
t8_dtri_face_child_face (const t8_dtri_t *triangle, int face, int face_child);

/** Given a face of an triangle return the face number of the parent of the triangle that matches the triangle's face.
 * Or return -1 if no face of the parent matches the face.
 * \param [in]  triangle The triangle.
 * \param [in]  face    Then number of the face.
 * \return              If \a face of \a elem is also a face of \a elem's parent, the face number of this face.
 *                      Otherwise -1.
 */
int
t8_dtri_face_parent_face (const t8_dtri_t *triangle, int face);

/** Given a triangle and a face of this triangle. If the face lies on the tree boundary, return the face number of the
 * tree face. If not the return value is arbitrary.
 * \param [in] element        The triangle.
 * \param [in] face     The index of a face of \a element.
 * \return The index of the tree face that \a face is a subface of, if \a face is on a tree boundary. Any arbitrary
 *         integer if \a element is not at a tree boundary.
 * \note For boundary triangles, this function is the inverse of \ref t8_dtri_root_face_to_face
 */
int
t8_dtri_tree_face (t8_dtri_t *element, int face);

/** Given a triangle and a face of the root triangle. If the triangle lies on the tree boundary, return the
 * corresponding face number of the triangle. If not the return value is arbitrary.
 * \param [in] element        The triangle.
 * \param [in] root_face     The index of a face of the root element.
 * \return The index of the face of \a element that is a subface of \a root_face, if \a element is on the tree boundary.
 *         Any arbitrary integer if \a element is not at a tree boundary.
 * \note For boundary triangles, this function is the inverse of \ref t8_dtri_tree_face
 */
int
t8_dtri_root_face_to_face (t8_dtri_t *element, int root_face);

/** Suppose we have two trees that share a common triangle f. Given a triangle e that is a subface of f in one of the
 * trees and given the orientation of the tree connection, construct the face triangle of the respective tree neighbor
 * that logically coincides with e but lies in the coordinate system of the neighbor tree.
 *  \param [in] trianglein     The face triangle.
 *  \param [in,out] triangle2 On return the face triangle \a triangle1 with respective to the coordinate system of the
 *                            other tree.
 *  \param [in] orientation The orientation of the tree-tree connection. \see t8_cmesh_set_join
 *  \param [in] sign      Depending on the topological orientation of the two tree faces, either 0
 *                        (both faces have opposite orientation) or 1 (both faces have the same top. orientattion).
 *                        \ref t8_eclass_face_orientation
 *  \param [in] is_smaller_face Flag to declare whether \a triangle1 belongs to the smaller face. A face f of tree T is
 *                        smaller than f' of T' if either the eclass of T is smaller or if the classes are equal and
 *                        f<f'. The orientation is defined in relation to the smaller face.
 * \note \a trianglein and \a triangle2 may point to the same element.
 */
void
t8_dtri_transform_face (const t8_dtri_t *trianglein, t8_dtri_t *triangle2, int orientation, int sign,
                        int is_smaller_face);

/** Test if a triangle lies inside of the root triangle, that is the triangle of level 0, anchor node (0,0) and type 0.
 *  \param [in]     element Input triangle.
 *  \return true    If \a element lies inside of the root triangle.
 */
int
t8_dtri_is_inside_root (t8_dtri_t *element);

/** Compute whether a given triangle shares a given face with its root tree.
 * \param [in] element        The input triangle.
 * \param [in] face     A face of \a element.
 * \return              True if \a face is a subface of the triangle's root element.
 */
int
t8_dtri_is_root_boundary (const t8_dtri_t *element, int face);

/** Test if two triangles have the same coordinates, type and level.
 * \return true if \a element1 describes the same triangle as \a element2.
 */
int
t8_dtri_is_equal (const t8_dtri_t *element1, const t8_dtri_t *element2);

/** Test if two triangles are siblings.
 * \param [in] element1 First triangle to be tested.
 * \param [in] element2 Second triangle to be tested.
 * \return true if \a element1 is equal to or a sibling of \a element2.
 */
int
t8_dtri_is_sibling (const t8_dtri_t *element1, const t8_dtri_t *element2);

/** Test if a triangle is the parent of another triangle.
 * \param [in] element triangle to be tested.
 * \param [in] child Possible child triangle.
 * \return true if \a element is the parent of \a child.
 */
int
t8_dtri_is_parent (const t8_dtri_t *element, const t8_dtri_t *child);

/** Test if a triangle is an ancestor of another triangle.
 * \param [in] element triangle to be tested.
 * \param [in] child Descendent triangle.
 * \return true if \a element is equal to or an ancestor of \a child.
 */
int
t8_dtri_is_ancestor (const t8_dtri_t *element, const t8_dtri_t *child);

/** Computes the linear position of a triangle in a uniform grid.
 * \param [in] element  triangle whose id will be computed.
 * \param [in] level level of uniform grid to be considered.
 * \return Returns the linear position of this triangle on a grid of level \a level.
 * \note This id is not the Morton index.
 */
t8_linearidx_t
t8_dtri_linear_id (const t8_dtri_t *element, int level);

/**
 * Same as init_linear_id, but we only consider the subtree. Used for computing the index of a tetrahedron lying in a
 * pyramid
 * \param [in, out] element   Existing triangle whose data will be filled
 * \param [in] id            Index to be considered
 * \param [in] start_level   The level of the root of the subtree
 * \param [in] end_level     Level of uniform grid to be considered
 * \param [in] parenttype    The type of the parent.
 */
void
t8_dtri_init_linear_id_with_level (t8_dtri_t *element, t8_linearidx_t id, const int start_level, const int end_level,
                                   t8_dtri_type_t parenttype);

/** Initialize a triangle as the triangle with a given global id in a uniform refinement of a given level. *
 * \param [in,out] element  Existing triangle whose data will be filled.
 * \param [in] id     Index to be considered.
 * \param [in] level  level of uniform grid to be considered.
 */
void
t8_dtri_init_linear_id (t8_dtri_t *element, t8_linearidx_t id, int level);

/** Initialize a triangle as the root triangle (type 0 at level 0)
 * \param [in,out] element Existing triangle whose data will be filled.
 */
void
t8_dtri_init_root (t8_dtri_t *element);

/** Computes the successor of a triangle in a uniform grid of level \a level.
 * \param [in] element  triangle whose id will be computed.
 * \param [in,out] s Existing triangle whose data will be filled with the
 *                data of t's successor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dtri_successor (const t8_dtri_t *element, t8_dtri_t *s, int level);

/** Compute the first descendant of a triangle at a given level. This is the descendant of
 * the triangle in a uniform maxlevel refinement that has the smaller id.
 * \param [in] element        Triangle whose descendant is computed.
 * \param [in] level    A given level. Must be grater or equal to \a element's level.
 * \param [out] s       Existing triangle whose data will be filled with the data
 *                      of t's first descendant.
 */
void
t8_dtri_first_descendant (const t8_dtri_t *element, t8_dtri_t *s, int level);

/** Compute the last descendant of a triangle at a given level. This is the descendant of
 * the triangle in a uniform maxlevel refinement that has the biggest id.
 * \param [in] element        Triangle whose descendant is computed.
 * \param [in] level    A given level. Must be grater or equal to \a element's level.
 * \param [out] s       Existing triangle whose data will be filled with the data
 *                      of t's last descendant.
 */
void
t8_dtri_last_descendant (const t8_dtri_t *element, t8_dtri_t *s, int level);

/** Compute the descendant of a triangle in a given corner.
 * \param [in] element       Triangle whose descendant is computed.
 * \param [out] s       Existing triangle whose data will be filled with the data
 *                      of t's descendant in \a corner.
 * \param [in]  corner  The corner in which the descendant should lie.
 * \param [in]  level   The refinement level of the descendant. Must be greater or
 *                      equal to \a element's level.
 */
void
t8_dtri_corner_descendant (const t8_dtri_t *element, t8_dtri_t *s, int corner, int level);

/** Computes the predecessor of a triangle in a uniform grid of level \a level.
 * \param [in] element  triangle whose id will be computed.
 * \param [in,out] s Existing triangle whose data will be filled with the
 *                data of \a element's predecessor on level \a level.
 * \param [in] level level of uniform grid to be considered.
 */
void
t8_dtri_predecessor (const t8_dtri_t *element, t8_dtri_t *s, int level);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] element  triangle to be considered.
 * \param [in] level level to be considered.
 * \return Returns its child id in 0..3
 */
int
t8_dtri_ancestor_id (const t8_dtri_t *element, int level);

/** Compute the position of the ancestor of this child at level \a level within its siblings.
 * \param [in] element  triangle to be considered.
 * \return Returns its child id in 0..3
 */
int
t8_dtri_child_id (const t8_dtri_t *element);

/** Return the level of a triangle.
 * \param [in] element  triangle to be considered.
 * \return        The level of \a element.
 */
int
t8_dtri_get_level (const t8_dtri_t *element);

/** Query whether all entries of a triangle are in valid ranges.
 * \param [in] element  triangle to be considered.
 * \return        True, if \a element is a valid triangle and it is safe to call any function on \a element. False otherwise.
 */
int
t8_dtri_is_valid (const t8_dtri_t *element);

#if T8_ENABLE_DEBUG
/** Set sensible default values for a triangle.
 * \param [in,out] element A triangle.
 */
void
t8_dtri_init (t8_dtri_t *element);
#endif

/** Packs an array of triangle elements into contiguous memory.
 * Triangles are packed as x and y coordinates, type and level.
 * Compare MPI_Pack function.
 * \param [in] elements The element array to be packed.
 * \param [in] count Number of elements to be packed.
 * \param [out] send_buffer Output buffer.
 * \param [in] buffer_size Output buffer size.
 * \param [in,out] position Current position in buffer.
 * \param [in] comm Communicator.
 */
void
t8_dtri_element_pack (t8_dtri_t **const elements, const unsigned int count, void *send_buffer, const int buffer_size,
                      int *position, sc_MPI_Comm comm);

/** Returns the upper bound on the amount of space needed to pack \a count elements.
 * Compare MPI_Pack_size function.
 * \param [in] count Number of packed elements.
 * \param [in] comm Communicator.
 * \param [out] pack_size Upper bound on size of packed elements.
 */
void
t8_dtri_element_pack_size (const unsigned int count, sc_MPI_Comm comm, int *pack_size);

/** Unpack a buffer of triangle elements into contiguous memory.
 * Compare MPI_Unpack function.
 * \param [in] recvbuf Buffer to be unpacked.
 * \param [in] buffer_size Size of \a recvbuf.
 * \param [in,out] position Current position in buffer.
 * \param [out] elements Array with the triangle elements (output buffer).
 * \param [in] count Number of elements to be unpacked.
 * \param [in] comm Communicator.
 */
void
t8_dtri_element_unpack (void *recvbuf, const int buffer_size, int *position, t8_dtri_t **elements,
                        const unsigned int count, sc_MPI_Comm comm);

T8_EXTERN_C_END ();

#endif /* T8_DTRI_BITS_H */

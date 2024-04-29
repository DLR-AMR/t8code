/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file p8est_bits.h
 *
 * Routines for manipulating quadrants (neighbors, parents, children, etc.)
 *
 * \ingroup p8est
 */

#ifndef P8EST_BITS_H
#define P8EST_BITS_H

#include <p8est.h>
#include <sc_random.h>

SC_EXTERN_C_BEGIN;

/** Write -1 into the pad8 and pad16 members of a quadrant.
 * This helps with valgrind cleanliness if a quadrant is sent over MPI.
 */
void                p8est_quadrant_pad (p8est_quadrant_t * q);

/** Prints one line with quadrant's x, y, z and level.
 * \param [in] log_priority  see \ref logpriorities in sc.h for the meanings
 *                           of numerical priority values
 * \param [in] q             quadrant to print
 */
void                p8est_quadrant_print (int log_priority,
                                          const p8est_quadrant_t * q);

/** Test if two quadrants have equal Morton indices.
 * \return true if \a q1 describes the same quadrant as \a q2.
 */
int                 p8est_quadrant_is_equal (const p8est_quadrant_t * q1,
                                             const p8est_quadrant_t * q2);

/** Copy the Morton indices of the quadrant \a q.
 *  \param[in] q      	 An extended quadrant.
 *  \param[in,out] copy  An existing quadrant that Morton indices will
 *                       be set to the Morton indices of \a q.
 */
void                p8est_quadrant_copy (const p8est_quadrant_t * q,
                                         p8est_quadrant_t * copy);

/** Test if two quadrants overlap.
 * \return true if \a q1 and \a q2 are equal or one is the ancestor of the
 * other.
 */
int                 p8est_quadrant_overlaps (const p8est_quadrant_t * q1,
                                             const p8est_quadrant_t * q2);

/** Test if two quadrants have equal Morton indices and the same tree id.
 * \return          true if \a q1 describes the same quadrant as \a q2
 *                  and the p.which_tree fields are equal.
 */
int                 p8est_quadrant_is_equal_piggy (const p8est_quadrant_t *
                                                   q1,
                                                   const p8est_quadrant_t *
                                                   q2);

/** Compare two quadrants in their Morton ordering.
 * Both quadrants must be valid.
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_quadrant_compare (const void *v1, const void *v2);

/** Compare two sets of coordinates in their Morton ordering.
 * Coordinates are signed, but the sorted order will treat them
 * as unsigned, with negative coordinates being greater than
 * positive coordinates because of their representation in twos-complement.
 * \param [in] v1, v2    Two sets of 3d coordinates.
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_coordinates_compare (const p4est_qcoord_t v1[],
                                               const p4est_qcoord_t v2[]);

/** Compare two quadrants in their Morton ordering, with equivalence if the
 * two quadrants overlap.
 * \return Returns < 0 if \a v1 < \a v2 and \a v1 and \a v2 do not overlap,
 *                   0 if \a v1 and \a v2 overlap,
 *                 > 0 if \a v1 > \a v2 and \a v1 and \a v2 do not overlap.
 */
int                 p8est_quadrant_disjoint (const void *v1, const void *v2);

/** Compare two quadrants in their Morton ordering and the which_tree member.
 * Both quadrants must be extended (superset of valid, see below).
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_quadrant_compare_piggy (const void *v1,
                                                  const void *v2);

/** Compare two quadrants with respect to their local_num in the piggy3 member.
 * \return Returns < 0 if \a v1 < \a v2,
 *                   0 if \a v1 == \a v2,
 *                 > 0 if \a v1 > \a v2
 */
int                 p8est_quadrant_compare_local_num (const void *v1,
                                                      const void *v2);

/** Test if two quadrants have equal Morton indices, callback version.
 * \return true if \a v1 describes the same quadrant as \a v2.
 */
int                 p8est_quadrant_equal_fn (const void *v1, const void *v2,
                                             const void *u);

/** Computes a hash value for a quadrant by the lookup3 method.
 */
unsigned            p8est_quadrant_hash_fn (const void *v, const void *u);

/** Test if two nodes are in the same tree and have equal Morton indices.
 * \param [in] v1   Pointer to a clamped or unclamped node, depending on u.
 * \param [in] v2   Pointer to a clamped or unclamped node, depending on u.
 * \param [in] u    User data, points to an int holding the clamped-flag.
 */
int                 p8est_node_equal_piggy_fn (const void *v1,
                                               const void *v2, const void *u);

/** Compute hash value of a node based on its tree and Morton index.
 * \param [in] v    Pointer to a clamped or unclamped node, depending on u.
 * \param [in] u    User data, points to an int holding the clamped-flag.
 */
unsigned            p8est_node_hash_piggy_fn (const void *v, const void *u);

/** Clamp a node inside the unit tree if it sits on a high border.
 * \param [in] n    Node to be clamped. Must not yet be clamped.
 * \param [out] r   Existing node overwritten by the clamped result.
 */
void                p8est_node_clamp_inside (const p8est_quadrant_t * n,
                                             p8est_quadrant_t * r);

/** Move a clamped node out on the border.
 * \param [in] n    Node to be unclamped in-place.
 */
void                p8est_node_unclamp (p8est_quadrant_t * n);

/** Find the enclosing quadrant of a given node at a given level.
 * \param [in] n        Clamped node.
 * \param [in] level    Level of the quadrant to be created.
 * \param [out] q       Output quadrant, n == q is permitted.
 */
void                p8est_node_to_quadrant (const p8est_quadrant_t * n,
                                            int level, p8est_quadrant_t * q);

/** Decide if a node is completely contained within a quadrant.
 * \param [in] q        Valid quadrant.
 * \param [in] n        Clamped node.
 */
int                 p8est_quadrant_contains_node (const p8est_quadrant_t * q,
                                                  const p8est_quadrant_t * n);

/** Compute the position of the ancestor of this child at level \a level within
 * its siblings.
 * \return Returns its child id in 0..7
 */
int                 p8est_quadrant_ancestor_id (const p8est_quadrant_t * q,
                                                int level);

/** Compute the position of this child within its siblings.
 * \return Returns its child id in 0..7
 */
int                 p8est_quadrant_child_id (const p8est_quadrant_t * q);

/** Test if Morton indices are inside the unit tree.
 * \param [in] coord   3d coordinates.
 * \return Returns true if \a (coord[0],coord[1],coord[2]) is inside the unit tree.
 */
int                 p8est_coordinates_is_inside_root (const p4est_qcoord_t
                                                      coord[]);

/** Test if a quadrant is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is inside the unit tree.
 */
int                 p8est_quadrant_is_inside_root (const p8est_quadrant_t *
                                                   q);

/** Test if a quadrant is inside the 3x3 box around the root tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is inside the unit tree.
 */
int                 p8est_quadrant_is_inside_3x3 (const p8est_quadrant_t * q);

/** Test if a quadrant is outside a tree face boundary (no edge or corner).
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree face.
 */
int                 p8est_quadrant_is_outside_face (const p8est_quadrant_t *
                                                    q);

/** Test if a quadrant is outside a tree edge boundary (no corner).
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree edge.
 */
int                 p8est_quadrant_is_outside_edge (const p8est_quadrant_t *
                                                    q);

/** Test if a quadrant is outside a tree edge boundary (no corner).
 * \param [in] q       Quadrant to be tested.
 * \param [out] edge   The tree edge number is computed if outside edge.
 *                     This pointer may be NULL.
 * \return Returns true if \a q is outside across a unit tree edge.
 */
int                 p8est_quadrant_is_outside_edge_extra (const
                                                          p8est_quadrant_t *
                                                          q, int *edge);

/** Test if a quadrant is outside a tree corner boundary.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is outside across a unit tree corner.
 */
int                 p8est_quadrant_is_outside_corner (const p8est_quadrant_t *
                                                      q);

/** Test if a quadrant is used to represent a mesh node.
 * \param [in] q        Quadrant to be tested.
 * \param [in] inside   If true, boundary nodes must be clamped inside.
 *                      If false, nodes must align with the quadrant grid.
 * \return Returns true if \a q is a node.
 */
int                 p8est_quadrant_is_node (const p8est_quadrant_t * q,
                                            int inside);

/** Test if Morton indices are valid and are inside the unit tree.
 * \param [in] coord  3d coordinates.
 * \param [in] level  level
 * \return Returns true if \a (coord[0],coord[1],coord[2],level) is valid.
 */
int                 p8est_coordinates_is_valid (const p4est_qcoord_t coord[],
                                                int level);

/** Test if a quadrant has valid Morton indices and is inside the unit tree.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is valid.
 */
int                 p8est_quadrant_is_valid (const p8est_quadrant_t * q);

/** Test if a quadrant has valid Morton indices in the 3x3 box around root.
 * \param [in] q Quadrant to be tested.
 * \return Returns true if \a q is extended.
 */
int                 p8est_quadrant_is_extended (const p8est_quadrant_t * q);

/** Test if two quadrants are siblings.
 * \param [in] q1 First quadrant to be tested.
 * \param [in] q2 Second quadrant to be tested.
 * \return true if \a q1 is unequal to and a sibling of \a q2.
 */
int                 p8est_quadrant_is_sibling (const p8est_quadrant_t * q1,
                                               const p8est_quadrant_t * q2);

/** Compute a specific child of a quadrant.
 * \param [in]     q    Input quadrant.
 * \param [in,out] r    Existing quadrant whose Morton index will be filled
 *                      with the coordinates of its child no. \b child_id.
 * \param [in] child_id The id of the child computed, 0..7.
 */
void                p8est_quadrant_child (const p8est_quadrant_t * q,
                                          p8est_quadrant_t * r, int child_id);

/** Test if two quadrants are siblings.
 * Descriptive, slower version of \a p8est_quadrant_is_sibling.
 * For debugging and educational purposes only.
 */
int                 p8est_quadrant_is_sibling_D (const p8est_quadrant_t * q1,
                                                 const p8est_quadrant_t * q2);

/** Test if 8 quadrants are siblings in Morton ordering.
 */
int                 p8est_quadrant_is_family (const p8est_quadrant_t * q0,
                                              const p8est_quadrant_t * q1,
                                              const p8est_quadrant_t * q2,
                                              const p8est_quadrant_t * q3,
                                              const p8est_quadrant_t * q4,
                                              const p8est_quadrant_t * q5,
                                              const p8est_quadrant_t * q6,
                                              const p8est_quadrant_t * q7);

/** Test if 8 quadrants are siblings in Morton ordering, array version.
 * \param [in] q   Array of 8 quadrants.
 */
int                 p8est_quadrant_is_familyv (const p8est_quadrant_t q[]);

/** Test if 8 quadrants are siblings in Morton ordering, array version.
 * \param [in] q   Array of 8 pointers to quadrants.
 */
int                 p8est_quadrant_is_familypv (p8est_quadrant_t * q[]);

/** Test if a quadrant is the parent of another quadrant.
 * \param [in] q Quadrant to be tested.
 * \param [in] r Possible child quadrant.
 * \return true if \a q is the parent of \a r.
 */
int                 p8est_quadrant_is_parent (const p8est_quadrant_t * q,
                                              const p8est_quadrant_t * r);

/** Test if a quadrant is the parent of another quadrant.
 * Descriptive, slower version of \a p8est_quadrant_is_parent.
 * For debugging and educational purposes only.
 */
int                 p8est_quadrant_is_parent_D (const p8est_quadrant_t * q,
                                                const p8est_quadrant_t * r);

/** Test if a quadrant is an ancestor of another quadrant.
 * \param [in] q Quadrant to be tested.
 * \param [in] r Descendent quadrant.
 * \return true if \a q is unequal to and an ancestor of \a r.
 */
int                 p8est_quadrant_is_ancestor (const p8est_quadrant_t * q,
                                                const p8est_quadrant_t * r);

/** Test if a quadrant is an ancestor of another quadrant.
 * Descriptive, slower version of \a p8est_quadrant_is_ancestor.
 * Contrary to \a p8est_quadrant_is_ancestor, it aborts for inter-tree q, r.
 * For debugging and educational purposes only.
 */
int                 p8est_quadrant_is_ancestor_D (const p8est_quadrant_t * q,
                                                  const p8est_quadrant_t * r);

/** Test if two quadrants follow each other in the tree with no holes.
 * \param [in] q A quadrant
 * \param [in] r Another quadrant
 * \return true if \a q is immediately before \a r in the tree.
 * \note for every \a q there are between 0 and P8EST_MAXLEVEL+1 possible nexts.
 */
int                 p8est_quadrant_is_next (const p8est_quadrant_t * q,
                                            const p8est_quadrant_t * r);

/** Test if two quadrants follow each other in the tree with no holes.
 * Descriptive, slower version of \a p8est_quadrant_is_next.
 * For debugging and educational purposes only.
 */
int                 p8est_quadrant_is_next_D (const p8est_quadrant_t * q,
                                              const p8est_quadrant_t * r);

/** Test if a quadrant has at least partial overlap with a tree.
 */
int                 p8est_quadrant_overlaps_tree (p8est_tree_t * tree,
                                                  const p8est_quadrant_t * q);

/** Test if a quadrant is completely contained within a tree.
 */
int                 p8est_quadrant_is_inside_tree (p8est_tree_t * tree,
                                                   const p8est_quadrant_t *
                                                   q);

/** Whether two descendants of a quadrant are first and last, up to size.
 * \param [in] f    An extended quadrant, need not be of maximum level.
 * \param [in] l    An extended quadrant, need not be of maximum level.
 *                  It must be greater equal \b f in the space filling curve.
 * \param [in] a    An extended quadrant,
 *                  equal to or ancestor of \b f, and likewise to/of \b l.
 * \return          Whether the first corner of \b f equals that of \b a and
 *                  the last corner of \b l equals that of \b a.
 */
int                 p8est_quadrant_is_first_last (const p8est_quadrant_t * f,
                                                  const p8est_quadrant_t * l,
                                                  const p8est_quadrant_t * a);

/** Enlarge a quadrant as long as its first corner stays the same.
 * We limit the enlargement by containing it in an ancestor quadrant.
 * \param [in] a    Extended quadrant.  On input and output, equal to or
 *                  strict ancestor of the quadrant \b q to be modified.
 * \param [in,out] q    On input and output, an extended quadrant and also
 *                      equal or a strict descendant of \b a.
 *                      Possibly enlarged by this function.
 */
void                p8est_quadrant_enlarge_first (const p8est_quadrant_t * a,
                                                  p8est_quadrant_t * q);

/** Enlarge a quadrant as long as its last corner stays the same.
 * We limit the enlargement by containing it in an ancestor quadrant.
 * \param [in] a    Extended quadrant.  On input and output, equal to or
 *                  strict ancestor of the quadrant \b q to be modified.
 * \param [in,out] q    On input and output, an extended quadrant and also
 *                      equal or a strict descendant of \b a.
 *                      Possibly enlarged by this function.
 */
void                p8est_quadrant_enlarge_last (const p8est_quadrant_t * a,
                                                 p8est_quadrant_t * q);

/** Compute the ancestor of a quadrant at a given level.
 * \param [in]  q       Input quadrant.
 * \param [in]  level   A smaller level than q.
 * \param [in,out]  r   Existing quadrant whose Morton index will be filled
 *                      with the ancestor of q at the given level.
 * \note The quadrant q may point to the same quadrant as r.
 *       The user_data of r are never modified.
 */
void                p8est_quadrant_ancestor (const p8est_quadrant_t * q,
                                             int level, p8est_quadrant_t * r);

/** Compute the parent of a quadrant.
 * \param [in]  q Input quadrant.
 * \param [in,out] r Existing quadrant whose Morton index will be filled
 *                   with the Morton index of the parent of \a q.
 *                   Its user_data will be untouched.
 * \note \a q may point to the same quadrant as \a r.
         The user_data of \a r is never modified.
 */
void                p8est_quadrant_parent (const p8est_quadrant_t * q,
                                           p8est_quadrant_t * r);

/** Compute a specific sibling of a quadrant.
 * \param [in]     q  Input quadrant.
 * \param [in,out] r  Existing quadrant whose Morton index will be filled
 *                    with the coordinates of sibling no. sibling_id of q.
 * \param [in]     sibling_id The id of the sibling computed, 0..3.
 */
void                p8est_quadrant_sibling (const p8est_quadrant_t * q,
                                            p8est_quadrant_t * r,
                                            int sibling_id);

/** Compute the face neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_face_neighbor (const p8est_quadrant_t * q,
                                                  int face,
                                                  p8est_quadrant_t * r);

/** Compute the face neighbor of a quadrant, transforming across tree
 * boundaries if necessary.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \a q.
 * \param [in]     face   The face across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 *                        By convention, if there is no tree across \a face,
 *                        \a r has the same Morton index as \a q.
 * \param [in,out] nface  if not NULL, set to the face of \a r that neighbors
 *                        \a q.  nface is encoded with orientation information
 *                        in the same manner as the tree_to_face array in
 *                        the p8est_connectivity_t struct.
 * \param [in]     conn   The connectivity structure for the forest.
 * \return Returns the tree that contains \a r.  By convention, if there is no
 * tree across \a face, then -1 is returned.
 */
p4est_locidx_t      p8est_quadrant_face_neighbor_extra (const p8est_quadrant_t
                                                        * q, p4est_topidx_t t,
                                                        int face,
                                                        p8est_quadrant_t * r,
                                                        int *nface,
                                                        p8est_connectivity_t *
                                                        conn);

/** Get the smaller face neighbors of \a q.
 *
 * Gets the smaller face neighbors, which are half of the size assuming the
 * 2-1 constant.
 *
 * The order of the \a n[i] is given in the Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the four smaller face neighbors.
 * \param [out] nur[0]..nur[3] If not NULL, filled with smallest quadrants
 *                     that fit in the upper right corners of \a n.
 */
void                p8est_quadrant_half_face_neighbors (const p8est_quadrant_t
                                                        * q, int face,
                                                        p8est_quadrant_t n[],
                                                        p8est_quadrant_t
                                                        nur[]);

/** Create all possible face neighbors of \a q.
 *
 * Gets the face neighbors, possible assuming the 2-1 constraint.
 * If the larger or smaller quadrants do not exist than they are returned
 * as initialized by P4EST_QUADRANT_INIT.
 *
 * The order of \a n[0] through \a n[3] are given in Morton ordering.
 *
 * \param [in]  q      The quadrant whose face neighbors will be constructed.
 * \param [in]  face   The face across which to generate the neighbors.
 * \param [out] n[0]..n[3] Filled with the smaller possible face neighbors,
 *                     which are half of the size if they exist
 *                     or initialized to P4EST_QUADRANT_INIT.
 * \param [out] n[4]   Filled with the face neighbor, which is the same size.
 * \param [out] n[5]   Filled with the face neighbor, which is twice the size
 *                     if it exists or initialized to P4EST_QUADRANT_INIT.
 */
void                p8est_quadrant_all_face_neighbors (const p8est_quadrant_t
                                                       * q, int face,
                                                       p8est_quadrant_t n[]);

/** Compute the edge neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     edge   The edge across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_edge_neighbor (const p8est_quadrant_t * q,
                                                  int edge,
                                                  p8est_quadrant_t * r);

/** Compute the edge neighbors of a quadrant, transforming across tree
 * boundaries if necessary.  Only computes neighbors that are not face
 * neighbors.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \a q.
 * \param [in]     edge   The edge across which to generate the neighbor.
 * \param [in,out] quads  An initialized but empty array where the edge
 *                        neighbors will be placed.
 * \param [in,out] treeids An initialized but empty array where the ids of the
 *                        trees containing the edge neighbors will be placed.
 * \param [in,out] nedges if not NULL, filled with the edges of \a quads that
 *                        neighbor \a q. the ints in \a nedges are encoded with
 *                        orientation information like the edge_to_edge array
 *                        in the p8est_connectivity_t struct
 * \param [in]     conn   The connectivity structure for the forest.
 */
void                p8est_quadrant_edge_neighbor_extra (const p8est_quadrant_t
                                                        * q, p4est_locidx_t t,
                                                        int edge, sc_array_t *
                                                        quads, sc_array_t *
                                                        treeids,
                                                        sc_array_t * nedges,
                                                        p8est_connectivity_t *
                                                        conn);

/** Compute the corner neighbor of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] r      Existing quadrant whose Morton index will be filled.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_corner_neighbor (const p8est_quadrant_t *
                                                    q, int corner,
                                                    p8est_quadrant_t * r);

/** Compute the corner neighbors of a quadrant, transforming across tree
 * boundaries if necessary.  Only computes neighbors that are not face or edge
 * neighbors.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     t      Tree that contains \a q.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] quads  An initialized but empty array where the corner
 *                        neighbors will be placed.
 * \param [in,out] treeids An initialized but empty array where the ids of the
 *                        trees containing the corner neighbors will be placed.
 * \param [in,out] ncorners if not NULL, filled with the corners of \a quads
 *                          that neighbor \a q.
 * \param [in]     conn   The connectivity structure for the forest.
 */
void                p8est_quadrant_corner_neighbor_extra (const
                                                          p8est_quadrant_t *
                                                          q, p4est_locidx_t t,
                                                          int corner,
                                                          sc_array_t * quads,
                                                          sc_array_t *
                                                          treeids,
                                                          sc_array_t *
                                                          ncorners,
                                                          p8est_connectivity_t
                                                          * conn);

/** Compute the half size corner neighbor of a quadrant.
 *
 * \param [in]  q       The quadrant whose corner neighbor will be constructed.
 * \param [in]  corner  The corner across which to generate the neighbor.
 * \param [out] r       Morton index filled with the half size corner neighbor.
 */
void                p8est_quadrant_half_corner_neighbor (const
                                                         p8est_quadrant_t * q,
                                                         int corner,
                                                         p8est_quadrant_t *
                                                         r);

/** Compute the corner node of a quadrant.
 * \param [in]     q      Input quadrant, must be valid.
 * \param [in]     corner The corner across which to generate the neighbor.
 * \param [in,out] r      Node that will not be clamped inside.
 * \note \a q may point to the same quadrant as \a r.
 */
void                p8est_quadrant_corner_node (const p8est_quadrant_t * q,
                                                int corner,
                                                p8est_quadrant_t * r);

/** Compute the 8 children of a quadrant.
 * \param [in]     q  Input quadrant.
 * \param [in,out] c0 First computed child.
 *                    \a q may point to the same quadrant as \a c0.
 * \note The user_data of \a c0, c1, c2, c3, c4, c5, c6, c7 is never modified.
 */
void                p8est_quadrant_children (const p8est_quadrant_t * q,
                                             p8est_quadrant_t * c0,
                                             p8est_quadrant_t * c1,
                                             p8est_quadrant_t * c2,
                                             p8est_quadrant_t * c3,
                                             p8est_quadrant_t * c4,
                                             p8est_quadrant_t * c5,
                                             p8est_quadrant_t * c6,
                                             p8est_quadrant_t * c7);

/** Compute the 8 children of a quadrant, array version.
 * \param [in]     q  Input quadrant.
 * \param [in,out] c  The 8 computed children in z-order.
 *                    q may point to the same quadrant as c[0].
 * \note The user_data of c[i] is never modified.
 */
void                p8est_quadrant_childrenv (const p8est_quadrant_t * q,
                                              p8est_quadrant_t c[]);

/** Compute the 8 children of a quadrant, array version.
 * \param [in]     q  Input quadrant.
 * \param [in,out] c  Pointers to the 8 computed children in z-order.
 *                    q may point to the same quadrant as c[0].
 * \note The user_data of c[i] is never modified.
 */
void                p8est_quadrant_childrenpv (const p8est_quadrant_t * q,
                                               p8est_quadrant_t * c[]);

/** Compute the first descendant of a quadrant on a given level.
 * \param [in]  q      Input quadrant.
 * \param [out] fd     First descendant of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void                p8est_quadrant_first_descendant (const p8est_quadrant_t *
                                                     q, p8est_quadrant_t * fd,
                                                     int level);

/** Compute the last descendant of a quadrant on a given level.
 * \param [in]  q      Input quadrant.
 * \param [out] ld     Last descendant of \a q on level \a level.
 * \param [in]  level  Level must be greater equal than q's level.
 */
void                p8est_quadrant_last_descendant (const p8est_quadrant_t *
                                                    q, p8est_quadrant_t * ld,
                                                    int level);

/** Compute the descendant of a quadrant touching a given corner.
 * \param [in]     q   Input quadrant.
 * \param [in,out] r   Existing quadrant whose Morton index will be filled.
 *                     Its user_data will be untouched.
 * \param [in]     c   The corner of \a q that \a r touches.
 * \param [in] level   The size of \a r.
 */
void                p8est_quadrant_corner_descendant (const p8est_quadrant_t *
                                                      q, p8est_quadrant_t * r,
                                                      int c, int level);

/** Computes the nearest common ancestor of two quadrants in the same tree.
 * \param [in]     q1 First input quadrant.
 * \param [in]     q2 Second input quadrant.
 * \param [in,out] r Existing quadrant whose Morton index will be filled.
 *                   Its user_data will be untouched.
 * \note \a q1, \a q2, \a r may point to the same quadrant.
 *       The user_data of \a r is never modified.
 */
void                p8est_nearest_common_ancestor (const p8est_quadrant_t *
                                                   q1,
                                                   const p8est_quadrant_t *
                                                   q2, p8est_quadrant_t * r);

/** Computes the nearest common ancestor of two quadrants in the same tree.
 * Descriptive, slower version of \a p8est_nearest_common_ancestor.
 * For debugging and educationial purposes only.
 */
void                p8est_nearest_common_ancestor_D (const p8est_quadrant_t *
                                                     q1,
                                                     const p8est_quadrant_t *
                                                     q2,
                                                     p8est_quadrant_t * r);

/** Transforms a quadrant/node across a face between trees.
 * \param [in]     q        Input quadrant/non-clamped node.
 * \param [in,out] r        Quadrant/node whose Morton index will be filled.
 * \param [in] ftransform   This array holds 9 integers.
 *             [0]..[2]     The coordinate axis sequence of the origin face.
 *             [3]..[5]     The coordinate axis sequence of the target face.
 *             [6]..[8]     Edge reverse flag for axes 0, 1; face code for 2.
 * \note \a q and \a r may NOT point to the same quadrant structure.
 */
void                p8est_quadrant_transform_face (const p8est_quadrant_t * q,
                                                   p8est_quadrant_t * r,
                                                   const int ftransform[]);

/** Transforms coordinates across a face between trees.
 * \param [in]  coords_in   Input coordinates.
 * \param [out] coords_out  Output coordinates.
 * \param [in] ftransform   This array holds 9 integers.
 *             [0]..[2]     The coordinate axis sequence of the origin face.
 *             [3]..[5]     The coordinate axis sequence of the target face.
 *             [6]..[8]     Edge reverse flag for axes 0, 1; face code for 2.
 */
void                p8est_coordinates_transform_face (const p4est_qcoord_t
                                                      coords_in[],
                                                      p4est_qcoord_t
                                                      coords_out[],
                                                      const int ftransform[]);

/** Checks if a quadrant touches an edge (diagonally inside or outside).
 */
int                 p8est_quadrant_touches_edge (const p8est_quadrant_t * q,
                                                 int edge, int inside);

/** Transforms a quadrant across an edge between trees.
 * \param [in]     q          Input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in]     ei         Edge info from p8est_find_edge_transform().
 * \param [in]     et         One of ei's transformations.
 * \param [in]     inside     The quadrant will be placed inside or outside.
 */
void                p8est_quadrant_transform_edge (const p8est_quadrant_t * q,
                                                   p8est_quadrant_t * r,
                                                   const p8est_edge_info_t *
                                                   ei,
                                                   const
                                                   p8est_edge_transform_t *
                                                   et, int inside);

/** Transforms coordinates on an edge between trees.
 * \param [in]     coords_in  Input coordinates.
 * \param [out]    coords_out Output coordinates.
 * \param [in]     ei         Edge info from p8est_find_edge_transform().
 * \param [in]     et         One of ei's transformations.
 */
void                p8est_coordinates_transform_edge (const p4est_qcoord_t
                                                      coords_in[],
                                                      p4est_qcoord_t
                                                      coords_out[],
                                                      const
                                                      p8est_edge_info_t *
                                                      ei,
                                                      const
                                                      p8est_edge_transform_t *
                                                      et);

/** Shifts a quadrant until it touches the specified edge from the inside.
 * If this shift is meant to recreate the effects of \a q on balancing across
 * the edge, then \a r, \a rup, and \a rdown may all be necessary for that
 * recreation.
 * \param [in]     q          Valid input quadrant.
 * \param [out]    r          Quadrant whose Morton index will be filled.
 *                            This quadrant results from shifting \a q
 *                            laterally towards the edge.
 * \param [out]    rup        Quadrant whose Morton index will be filled (may
 *                            be NULL).  This quadrant results from shifting
 *                            \a q diagonally towards \a edge's higher
 *                            corner.
 * \param [out]    rdown      Quadrant whose Morton index will be filled (may
 *                            be NULL).  This quadrant results from shifting
 *                            \a q diagonally towards \a edge's lower corner.
 * \param [in]     edge       Edge index.
 */
void                p8est_quadrant_shift_edge (const p8est_quadrant_t * q,
                                               p8est_quadrant_t * r,
                                               p8est_quadrant_t * rup,
                                               p8est_quadrant_t * rdown,
                                               int edge);

/** Checks if a quadrant touches a corner (diagonally inside or outside).
 */
int                 p8est_quadrant_touches_corner (const p8est_quadrant_t * q,
                                                   int corner, int inside);

/** Move a quadrant inside or diagonally outside a corner position.
 * \param [in,out] q        This quadrant only requires a valid level.
 * \param [in]     icorner  Number of the corner in 0..7.
 * \param [int]    inside   Boolean flag for inside or diagonally outside.
 */
void                p8est_quadrant_transform_corner (p8est_quadrant_t * r,
                                                     int corner, int inside);

/** Shifts a quadrant until it touches the specified corner from the inside.
 * \param [in]     q          Valid input quadrant.
 * \param [in,out] r          Quadrant whose Morton index will be filled.
 * \param [in]     corner     Corner index.
 */
void                p8est_quadrant_shift_corner (const p8est_quadrant_t * q,
                                                 p8est_quadrant_t * r,
                                                 int corner);

/** Computes the linear position of a quadrant in a uniform grid.
 * The grid and quadrant levels need not coincide.
 * If they do, this is the inverse of \ref p8est_quadrant_set_morton.
 * \param [in] quadrant  Quadrant whose linear index will be computed.
 *                       If the quadrant is smaller than the grid (has a higher
 *                       quadrant->level), the result is computed from its
 *                       ancestor at the grid's level.
 *                       If the quadrant has a smaller level than the grid (it
 *                       is bigger than a grid cell), the grid cell sharing its
 *                       lower left corner is used as reference.
 * \param [in] level     The level of the regular grid compared to which the
 *                       linear position is to be computed.
 *                       The level must be less equal P8EST_OLD_MAXLEVEL
 *                       since this is a legacy function restricted to 64 bits.
 * \return Returns the linear position of this quadrant on a grid.
 * \note The user_data of \a quadrant is never modified.
 */
uint64_t            p8est_quadrant_linear_id (const p8est_quadrant_t *
                                              quadrant, int level);

/** Set quadrant Morton indices based on linear position in uniform grid.
 * This is the inverse operation of \ref p8est_quadrant_linear_id.
 * \param [in,out] quadrant  Quadrant whose Morton indices will be set.
 * \param [in]     level     Level of the grid and of the resulting quadrant.
 * \param [in]     id        Linear index of the quadrant on a uniform grid.
 * \note The user_data of \a quadrant is never modified.
 */
void                p8est_quadrant_set_morton (p8est_quadrant_t * quadrant,
                                               int level, uint64_t id);

/** Compute the successor according to the Morton index in a uniform mesh.
 * \param[in] quadrant  Quadrant whose Morton successor will be computed.
 *                      Must not be the last (top right) quadrant in the tree.
 * \param[in,out] result    The coordinates and level of the successor of
 *                          \b quadrant will be saved in \b result.
 */
void                p8est_quadrant_successor (const p8est_quadrant_t *
                                              quadrant,
                                              p8est_quadrant_t * result);

/** Compute the predecessor according to the Morton index in a uniform mesh.
 * \param[in] quadrant  Quadrant whose Morton predecessor will be computed.
 *                      Must not be the first (bottom left) quadrant in the tree.
 * \param[in,out] result    The coordinates and level of the predecessor of
 *                          \b quadrant will be saved in result.
 */
void                p8est_quadrant_predecessor (const p8est_quadrant_t *
                                                quadrant,
                                                p8est_quadrant_t * result);

/** Initialize a random number generator by quadrant coordinates.
 * This serves to generate partition-independent and reproducible samples.
 * \param [in]  q               Valid quadrant.
 * \param [out] rstate          New state of random number generator.
 */
void                p8est_quadrant_srand (const p8est_quadrant_t * q,
                                          sc_rand_state_t * rstate);

/** Transform a quadrant from self's coordinate system to neighbor's coordinate system.
 *
 * \param [in]  nt            A neighbor transform.
 * \param [in]  self_quad     Input quadrant in self coordinates.
 * \param [out] neigh_coords  Quad transformed into neighbor coordinates.
 *
 * \note This transform gives meaningful results when \a self_quad is inside
 * the tree root or touches the interface between the two trees in the
 * transform.
 */
void                p8est_neighbor_transform_quadrant
  (const p8est_neighbor_transform_t * nt,
   const p8est_quadrant_t * self_quad, p8est_quadrant_t * neigh_quad);

/** Transform a quadrant from a neighbors's coordinate system to self's coordinate system.
 *
 * \param [in]  nt            A neighbor transform.
 * \param [in]  neigh_coords  Input quadrant in neighbor coordinates.
 * \param [out] self_coords   Quad transformed into self coordinates.
 *
 * \note This transform gives meaningful results when \a neigh_quad is inside
 * the tree root or touches the interface between the two trees in the
 * transform.
 */
void                p8est_neighbor_transform_quadrant_reverse
  (const p8est_neighbor_transform_t * nt,
   const p8est_quadrant_t * neigh_quad, p8est_quadrant_t * self_quad);

/** Check if a descendant shares a face with a (strict) ancestor.
 *
 * \param [in]  descendant   The descendant in question.
 * \param [in]  ancestor     The ancestor must not be equal to the descendant.
 * \param [in]  face         The face of the descendant.
 *
 * \return true if descendant face touches ancestor face, false otherwise.
*/
int                 p8est_quadrant_is_ancestor_face (const p8est_quadrant_t *
                                                     descendant,
                                                     const p8est_quadrant_t *
                                                     ancestor, int face);

/** Check if a descendant shares a corner with a (strict) ancestor.
 *
 * \param [in]  descendant   The descendant in question.
 * \param [in]  ancestor     The ancestor must not be equal to the descendant.
 * \param [in]  corner       The corner of the descendant.
 *
 * \return  true if descendant face touches ancestor corner, false otherwise.
*/
int                 p8est_quadrant_is_ancestor_corner (const p8est_quadrant_t
                                                       * descendant,
                                                       const p8est_quadrant_t
                                                       * ancestor,
                                                       int corner);

SC_EXTERN_C_END;

#endif /* !P8EST_BITS_H */

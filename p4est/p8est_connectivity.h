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

/** \file p8est_connectivity.h
 *
 * The connectivity defines the coarse topology of the forest.
 *
 * A 3D forest consists of one or more octrees, each of which a logical
 * cube.
 * Each tree has a local coordinate system, which defines the origin and the
 * direction of its x-, y-, and z-axes as well as the numbering of its
 * faces, edges, and corners.
 * Each tree may connect to any other tree (including itself) across any of
 * its faces, edges and/or corners, where the neighbor may be arbitrarily
 * rotated and/or flipped.
 * The \ref p8est_connectivity data structure stores these connections.
 *
 * We impose the following requirement for consistency of \ref p8est_balance :
 *
 * \note If a connectivity implies natural connections between trees that
 * are edge neighbors without being face neighbors, these edges shall be
 * encoded explicitly in the connectivity.  If a connectivity implies
 * natural connections between trees that are corner neighbors without being
 * edge or face neighbors, these corners shall be encoded explicitly in the
 * connectivity.
 * Please see the documentation of \ref p8est_connectivity_t for the exact
 * encoding convention.
 *
 * We provide various predefined connectivities by dedicated constructors,
 * such as
 *
 *  * \ref p8est_connectivity_new_unitcube for the unit square,
 *  * \ref p8est_connectivity_new_periodic for the periodic unit square,
 *  * \ref p8est_connectivity_new_brick for a rectangular grid of trees,
 *  * \ref p8est_connectivity_new_drop for a sparsely loop of trees.
 *
 * \ingroup p8est
 */

#ifndef P8EST_CONNECTIVITY_H
#define P8EST_CONNECTIVITY_H

#include <sc_io.h>
#include <p4est_base.h>

SC_EXTERN_C_BEGIN;

/** The spatial dimension */
#define P8EST_DIM 3
/** The number of faces of an octant
 *
 * \note for uniform naming reasons, an
 * octant is represented by the datatype p8est_quadrant_t */
#define P8EST_FACES (2 * P8EST_DIM)
/** The number of children of an octant
 *
 * also the nmber of corners */
#define P8EST_CHILDREN 8
/** The number of children/corners touching one face */
#define P8EST_HALF (P8EST_CHILDREN / 2)
/** The number of edges around an octant */
#define P8EST_EDGES 12
/** The size of insulation layer */
#define P8EST_INSUL 27

/** Only use logical AND term in 3D */
#define P8EST_ONLY_P8_LAND(x) && (x)

/** Only use comma and expression in 3D */
#define P8EST_ONLY_P8_COMMA(x) , (x)

/** Exponentiate with dimension */
#define P8EST_DIM_POW(a) ((a) * (a) * (a))

/** size of face transformation encoding */
#define P8EST_FTRANSFORM 9

/** p8est identification string */
#define P8EST_STRING "p8est"

/** Increase this number whenever the on-disk format for
 * p8est_connectivity, p8est, or any other 3D data structure changes.
 * The format for reading and writing must be the same.
 */
#define P8EST_ONDISK_FORMAT 0x3000009

/** Characterize a type of adjacency.
 *
 * Several functions involve relationships between neighboring trees and/or
 * quadrants, and their behavior depends on how one defines adjacency:
 * 1) entities are adjacent if they share a face, or
 * 2) entities are adjacent if they share a face or corner, or
 * 3) entities are adjacent if they share a face, corner or edge.
 * p8est_connect_type_t is used to choose the desired behavior.
 * This enum must fit into an int8_t.
 */
typedef enum
{
  /* make sure to have different values 2D and 3D */
  P8EST_CONNECT_SELF = 30,      /**< No balance whatsoever. */
  P8EST_CONNECT_FACE = 31,      /**< Balance across faces only. */
  P8EST_CONNECT_EDGE = 32,      /**< Balance across faces and edges. */
  P8EST_CONNECT_ALMOST = P8EST_CONNECT_EDGE,    /**< = CORNER - 1. */
  P8EST_CONNECT_CORNER = 33,    /**< Balance faces, edges, corners. */
  P8EST_CONNECT_FULL = P8EST_CONNECT_CORNER     /**< = CORNER. */
}
p8est_connect_type_t;

#ifdef P4EST_BACKWARD_DEALII
typedef p8est_connect_type_t p8est_balance_type_t;
#endif

/** Typedef for serialization method. */
typedef enum
{
  P8EST_CONN_ENCODE_NONE = SC_IO_ENCODE_NONE,
  P8EST_CONN_ENCODE_LAST        /**< Invalid entry to close the list. */
}
p8est_connectivity_encode_t;

/** Convert the p8est_connect_type_t into a number.
 * \param [in] btype    The balance type to convert.
 * \return              Returns 1, 2 or 3.
 */
int                 p8est_connect_type_int (p8est_connect_type_t btype);

/** Convert the p8est_connect_type_t into a const string.
 * \param [in] btype    The balance type to convert.
 * \return              Returns a pointer to a constant string.
 */
const char         *p8est_connect_type_string (p8est_connect_type_t btype);

/** This structure holds the 3D inter-tree connectivity information.
 * Identification of arbitrary faces, edges and corners is possible.
 *
 * The arrays tree_to_* are stored in z ordering.
 * For corners the order wrt. zyx is 000 001 010 011 100 101 110 111.
 * For faces the order is -x +x -y +y -z +z.
 * They are allocated [0][0]..[0][N-1]..[num_trees-1][0]..[num_trees-1][N-1].
 * where N is 6 for tree and face, 8 for corner, 12 for edge.
 * If a face is on the physical boundary it must connect to itself.
 *
 * The values for tree_to_face are in 0..23
 * where ttf % 6 gives the face number and ttf / 6 the face orientation code.
 * The orientation is determined as follows.  Let my_face and other_face
 * be the two face numbers of the connecting trees in 0..5.  Then the first
 * face corner of the lower of my_face and other_face connects to a face
 * corner numbered 0..3 in the higher of my_face and other_face.  The face
 * orientation is defined as this number.  If my_face == other_face, treating
 * either of both faces as the lower one leads to the same result.
 *
 * It is valid to specify num_vertices as 0.
 * In this case vertices and tree_to_vertex are set to NULL.
 * Otherwise the vertex coordinates are stored in the array vertices as
 * [0][0]..[0][2]..[num_vertices-1][0]..[num_vertices-1][2].
 * Vertex coordinates are optional and not used for inferring topology.
 *
 * The edges are stored when they connect trees that are not already face
 * neighbors at that specific edge.
 * In this case tree_to_edge indexes into \a ett_offset.
 * Otherwise the tree_to_edge entry must be -1 and this edge is ignored.
 * If num_edges == 0, tree_to_edge and edge_to_* arrays are set to NULL.
 *
 * The arrays edge_to_* store a variable number of entries per edge.
 * For edge e these are at position [ett_offset[e]]..[ett_offset[e+1]-1].
 * Their number for edge e is ett_offset[e+1] - ett_offset[e].
 * The entries encode all trees adjacent to edge e.
 * The size of the edge_to_* arrays is num_ett = ett_offset[num_edges].
 * The edge_to_edge array holds values in 0..23, where the lower 12 indicate
 * one edge orientation and the higher 12 the opposite edge orientation.
 *
 * The corners are stored when they connect trees that are not already edge
 * or face neighbors at that specific corner.
 * In this case tree_to_corner indexes into \a ctt_offset.
 * Otherwise the tree_to_corner entry must be -1 and this corner is ignored.
 * If num_corners == 0, tree_to_corner and corner_to_* arrays are set to NULL.
 *
 * The arrays corner_to_* store a variable number of entries per corner.
 * For corner c these are at position [ctt_offset[c]]..[ctt_offset[c+1]-1].
 * Their number for corner c is ctt_offset[c+1] - ctt_offset[c].
 * The entries encode all trees adjacent to corner c.
 * The size of the corner_to_* arrays is num_ctt = ctt_offset[num_corners].
 *
 * The *_to_attr arrays may have arbitrary contents defined by the user.
 *
 * \note If a connectivity implies natural connections between trees that
 * are edge neighbors without being face neighbors, these edges shall be
 * encoded explicitly in the connectivity.  If a connectivity implies
 * natural connections between trees that are corner neighbors without being
 * edge or face neighbors, these corners shall be encoded explicitly in the
 * connectivity.
 */
typedef struct p8est_connectivity
{
  p4est_topidx_t      num_vertices; /**< the number of vertices that define
                                         the \a embedding of the forest (not
                                         the topology) */
  p4est_topidx_t      num_trees;    /**< the number of trees */
  p4est_topidx_t      num_edges;    /**< the number of edges that help define
                                         the topology */
  p4est_topidx_t      num_corners;  /**< the number of corners that help
                                         define the topology */

  double             *vertices;     /**< an array of size
                                         (3 * \a num_vertices) */
  p4est_topidx_t     *tree_to_vertex; /**< embed each tree into \f$R^3\f$ for
                                           e.g. visualization (see
                                           p8est_vtk.h) */

  size_t              tree_attr_bytes;  /**< bytes per tree in tree_to_attr */
  char               *tree_to_attr;     /**< not touched by p4est */

  p4est_topidx_t     *tree_to_tree; /**< (6 * \a num_trees) neighbors across
                                         faces */
  int8_t             *tree_to_face; /**< (6 * \a num_trees) face to
                                         face+orientation (see description) */
  p4est_topidx_t     *tree_to_edge; /**< (12 * \a num_trees) or NULL (see
                                          description) */
  p4est_topidx_t     *ett_offset; /**< edge to offset in \a edge_to_tree and
                                       \a edge_to_edge */
  p4est_topidx_t     *edge_to_tree; /**< list of trees that meet at an edge */
  int8_t             *edge_to_edge; /**< list of tree-edges+orientations that
                                         meet at an edge (see description) */
  p4est_topidx_t     *tree_to_corner; /**< (8 * \a num_trees) or NULL (see
                                           description) */
  p4est_topidx_t     *ctt_offset; /**< corner to offset in \a corner_to_tree
                                       and \a corner_to_corner */
  p4est_topidx_t     *corner_to_tree; /**< list of trees that meet at a corner */
  int8_t             *corner_to_corner; /**< list of tree-corners that meet at
                                             a corner */
}
p8est_connectivity_t;

/** Calculate memory usage of a connectivity structure.
 * \param [in] conn   Connectivity structure.
 * \return            Memory used in bytes.
 */
size_t              p8est_connectivity_memory_used (p8est_connectivity_t *
                                                    conn);

/** Generic interface for transformations between a tree and any of its edge */
typedef struct
{
  p4est_topidx_t      ntree; /**< The number of the tree*/
  int8_t              nedge; /**< The number of the edge*/
  int8_t              naxis[3]; /**< The 3 edge coordinate axes*/
  int8_t              nflip; /**< The orientation of the edge*/
  int8_t              corners; /**< The corners connected to the edge*/
}
p8est_edge_transform_t;

/** Information about the neighbors of an edge*/
typedef struct
{
  int8_t              iedge; /**< The information of the edge*/
  sc_array_t          edge_transforms; /**< The array of neighbors of the originating
                                         edge */
}
p8est_edge_info_t;

/** Generic interface for transformations between a tree and any of its corner */
typedef struct
{
  p4est_topidx_t      ntree; /**< The number of the tree*/
  int8_t              ncorner; /**< The number of the corner*/
}
p8est_corner_transform_t;

/** Information about the neighbors of a corner*/
typedef struct
{
  p4est_topidx_t      icorner; /**< The number of the originating corner */
  sc_array_t          corner_transforms; /**< The array of neighbors of the originating
                                         corner */
}
p8est_corner_info_t;

/** Generic interface for transformations between a tree and any of its neighbors */
typedef struct
{
  p8est_connect_type_t neighbor_type; /**< type of connection to neighbor*/
  p4est_topidx_t      neighbor;   /**< neighbor tree index */
  int8_t              index_self; /**< index of interface from self's
                                       perspective */
  int8_t              index_neighbor; /**< index of interface from neighbor's
                                           perspective */
  int8_t              perm[P8EST_DIM]; /**< permutation of dimensions when
                                            transforming self coords to
                                            neighbor coords */
  int8_t              sign[P8EST_DIM]; /**< sign changes when transforming self
                                            coords to neighbor coords */
  p4est_qcoord_t      origin_self[P8EST_DIM]; /**< point on the interface from
                                                  self's perspective */
  p4est_qcoord_t      origin_neighbor[P8EST_DIM]; /**< point on the interface
                                                      from neighbor's
                                                      perspective */
}
p8est_neighbor_transform_t;

/* *INDENT-OFF* */

/** Transform from self's coordinate system to neighbor's coordinate system.
 *
 * \param [in]  nt            A neighbor transform.
 * \param [in]  self_coords   Input quadrant coordinates in self coordinates.
 * \param [out] neigh_coords  Coordinates transformed into neighbor coordinates.
 */
void                p8est_neighbor_transform_coordinates
                      (const p8est_neighbor_transform_t * nt,
                       const p4est_qcoord_t self_coords[P8EST_DIM],
                       p4est_qcoord_t neigh_coords[P8EST_DIM]);

/** Transform from neighbor's coordinate system to self's coordinate system.
 *
 * \param [in]  nt            A neighbor transform.
 * \param [in]  neigh_coords  Input quadrant coordinates in self coordinates.
 * \param [out] self_coords   Coordinates transformed into neighbor coordinates.
 */
void                p8est_neighbor_transform_coordinates_reverse
                      (const p8est_neighbor_transform_t * nt,
                       const p4est_qcoord_t neigh_coords[P8EST_DIM],
                       p4est_qcoord_t self_coords[P8EST_DIM]);

/**  Fill an array with the neighbor transforms based on a specific boundary type.
 *   This function generalizes all other inter-tree transformation objects
 *
 * \param [in]  conn   Connectivity structure.
 * \param [in]  tree_id The number of the tree.
 * \param [in]  boundary_type  Type of boundary connection (self, face, edge, corner).
 * \param [in]  boundary_index  The index of the boundary.
 * \param [in,out] neighbor_transform_array   Array of the neighbor transforms.
 */
void                p8est_connectivity_get_neighbor_transforms
                      (p8est_connectivity_t *conn,
                       p4est_topidx_t tree_id,
                       p8est_connect_type_t boundary_type,
                       int boundary_index,
                       sc_array_t *neighbor_transform_array);

/* *INDENT-ON* */

/** Store the corner numbers 0..7 for each tree face. */
extern const int    p8est_face_corners[6][4];

/** Store the edge numbers 0..12 for each tree face. */
extern const int    p8est_face_edges[6][4];

/** Store the face numbers in the face neighbor's system. */
extern const int    p8est_face_dual[6];

/* face corners */

/** Store only the 8 out of 24 possible permutations that occur. */
extern const int    p8est_face_permutations[8][4];

/** Store the 3 occurring sets of 4 permutations per face. */
extern const int    p8est_face_permutation_sets[3][4];

/** For each face combination store the permutation set.
 * The order is [my_face][neighbor_face] */
extern const int    p8est_face_permutation_refs[6][6];

/* face edges */

/** Store only the 8 out of 24 possible permutations that occur. */
extern const int    p8est_face_edge_permutations[8][4];

/** Store the 3 occurring sets of 4 permutations per face. */
extern const int    p8est_face_edge_permutation_sets[3][4];

/** Store the face numbers 0..5 adjacent to each tree edge. */
extern const int    p8est_edge_faces[12][2];

/** Store the corner numbers 0..8 for each tree edge. */
extern const int    p8est_edge_corners[12][2];

/** Store the edge corner numbers 0..1 for the corners touching a tree edge
    or -1 if combination is invalid */
extern const int    p8est_edge_edge_corners[12][8];

/** Store the face corner numbers 0..3 for the faces touching a tree edge.
    Is -1 for invalid combinations of indices */
extern const int    p8est_edge_face_corners[12][6][2];

/** Store the face edge numbers 0..3 for the faces touching a tree edge.
    Is -1 for invalid combinations of indices */
extern const int    p8est_edge_face_edges[12][6];

/** Store the face numbers 0..5 for each tree corner. */
extern const int    p8est_corner_faces[8][3];

/** Store the edge numbers 0..11 for each tree corner. */
extern const int    p8est_corner_edges[8][3];

/** Store the face corner numbers for the faces touching a tree corner.
    Is -1 for invalid combinations. */
extern const int    p8est_corner_face_corners[8][6];

/** Store the edge corner numbers for the edges touching a tree corner.
    Is -1 for invalid combinations. */
extern const int    p8est_corner_edge_corners[8][12];

/** Store the faces for each child and edge, can be -1. */
extern const int    p8est_child_edge_faces[8][12];

/** Store the faces for each child and corner, can be -1. */
extern const int    p8est_child_corner_faces[8][8];

/** Store the edges for each child and corner, can be -1. */
extern const int    p8est_child_corner_edges[8][8];

/** Transform a corner across one of the adjacent faces into a neighbor tree.
 * It expects a face permutation index that has been precomputed.
 * \param [in] c    A corner number in 0..7.
 * \param [in] f    A face number that touches the corner \a c.
 * \param [in] nf   A neighbor face that is on the other side of \a f.
 * \param [in] set  A value from \a p8est_face_permutation_sets that is
 *                  obtained using \a f, \a nf, and a valid orientation:
 *                  ref = p8est_face_permutation_refs[f][nf];
 *                  set = p8est_face_permutation_sets[ref][orientation];
 * \return          The corner number in 0..7 seen from the other face.
 */
int                 p8est_connectivity_face_neighbor_corner_set
  (int c, int f, int nf, int set);

/** Transform a face corner across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] fc   A face corner number in 0..3.
 * \param [in] f    A face that the face corner \a fc is relative to.
 * \param [in] nf   A neighbor face that is on the other side of \a f.
 * \param [in] o    The orientation between tree boundary faces \a f and \a nf.
 * \return          The face corner number relative to the neighbor's face.
 */
int                 p8est_connectivity_face_neighbor_face_corner
  (int fc, int f, int nf, int o);

/** Transform a corner across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] c    A corner number in 0..7.
 * \param [in] f    A face number that touches the corner \a c.
 * \param [in] nf   A neighbor face that is on the other side of \a f.
 * \param [in] o    The orientation between tree boundary faces \a f and \a nf.
 * \return          The number of the corner seen from the neighbor tree.
 */
int                 p8est_connectivity_face_neighbor_corner
  (int c, int f, int nf, int o);

/** Transform a face-edge across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] fe   A face edge number in 0..3.
 * \param [in] f    A face number that touches the edge \a e.
 * \param [in] nf   A neighbor face that is on the other side of \a f.
 * \param [in] o    The orientation between tree boundary faces \a f and \a nf.
 * \return          The face edge number seen from the neighbor tree.
 */
int                 p8est_connectivity_face_neighbor_face_edge
  (int fe, int f, int nf, int o);

/** Transform an edge across one of the adjacent faces into a neighbor tree.
 * This version expects the neighbor face and orientation separately.
 * \param [in] e    A edge number in 0..11.
 * \param [in] f    A face 0..5 that touches the edge \a e.
 * \param [in] nf   A neighbor face that is on the other side of \a f.
 * \param [in] o    The orientation between tree boundary faces \a f and \a nf.
 * \return          The edge's number seen from the neighbor.
 */
int                 p8est_connectivity_face_neighbor_edge
  (int e, int f, int nf, int o);

/** Transform an edge corner across one of the adjacent edges into a neighbor tree.
 * \param [in] ec   An edge corner number in 0..1.
 * \param [in] o    The orientation of a tree boundary edge connection.
 * \return          The edge corner number seen from the other tree.
 */
int                 p8est_connectivity_edge_neighbor_edge_corner
  (int ec, int o);

/** Transform a corner across one of the adjacent edges into a neighbor tree.
 * This version expects the neighbor edge and orientation separately.
 * \param [in] c    A corner number in 0..7.
 * \param [in] e    An edge 0..11 that touches the corner \a c.
 * \param [in] ne   A neighbor edge that is on the other side of \a e.
 * \param [in] o    The orientation between tree boundary edges \a e and \a ne.
 * \return          Corner number seen from the neighbor.
 */
int                 p8est_connectivity_edge_neighbor_corner
  (int c, int e, int ne, int o);

/** Allocate a connectivity structure.
 * The attribute fields are initialized to NULL.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_edges      Number of tree-connecting edges.
 * \param [in] num_ett        Number of total trees in edge_to_tree array.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] num_ctt        Number of total trees in corner_to_tree array.
 * \return                    A connectivity structure with allocated arrays.
 */
p8est_connectivity_t *p8est_connectivity_new (p4est_topidx_t num_vertices,
                                              p4est_topidx_t num_trees,
                                              p4est_topidx_t num_edges,
                                              p4est_topidx_t num_ett,
                                              p4est_topidx_t num_corners,
                                              p4est_topidx_t num_ctt);

/** Allocate a connectivity structure and populate from constants.
 * The attribute fields are initialized to NULL.
 * \param [in] num_vertices   Number of total vertices (i.e. geometric points).
 * \param [in] num_trees      Number of trees in the forest.
 * \param [in] num_edges      Number of tree-connecting edges.
 * \param [in] num_corners    Number of tree-connecting corners.
 * \param [in] vertices       Coordinates of the vertices of the trees.
 * \param [in] ttv            The tree-to-vertex array.
 * \param [in] ttt            The tree-to-tree array.
 * \param [in] ttf            The tree-to-face array (int8_t).
 * \param [in] tte            The tree-to-edge array.
 * \param [in] eoff           Edge-to-tree offsets (num_edges + 1 values).
 *                            This must always be non-NULL; in trivial cases
 *                            it is just a pointer to a p4est_topix value of 0.
 * \param [in] ett            The edge-to-tree array.
 * \param [in] ete            The edge-to-edge array.
 * \param [in] ttc            The tree-to-corner array.
 * \param [in] coff           Corner-to-tree offsets (num_corners + 1 values).
 *                            This must always be non-NULL; in trivial cases
 *                            it is just a pointer to a p4est_topix value of 0.
 * \param [in] ctt            The corner-to-tree array.
 * \param [in] ctc            The corner-to-corner array.
 * \return                    The connectivity is checked for validity.
 */
p8est_connectivity_t *p8est_connectivity_new_copy (p4est_topidx_t
                                                   num_vertices,
                                                   p4est_topidx_t num_trees,
                                                   p4est_topidx_t num_edges,
                                                   p4est_topidx_t num_corners,
                                                   const double *vertices,
                                                   const p4est_topidx_t * ttv,
                                                   const p4est_topidx_t * ttt,
                                                   const int8_t * ttf,
                                                   const p4est_topidx_t * tte,
                                                   const p4est_topidx_t *
                                                   eoff,
                                                   const p4est_topidx_t * ett,
                                                   const int8_t * ete,
                                                   const p4est_topidx_t * ttc,
                                                   const p4est_topidx_t *
                                                   coff,
                                                   const p4est_topidx_t * ctt,
                                                   const int8_t * ctc);

/** Broadcast a connectivity structure that exists only on one process to all.
 *  On the other processors, it will be allocated using p8est_connectivity_new.
 *  \param [in] conn_in For the root process the connectivity to be broadcast,
 *                      for the other processes it must be NULL.
 *  \param [in] root    The rank of the process that provides the connectivity.
 *  \param [in] comm    The MPI communicator.
 *  \return             For the root process this is a pointer to \a conn_in.
 *                      Else, a pointer to a newly allocated connectivity
 *                      structure with the same values as \a conn_in on the
 *                      root process.
 */
p8est_connectivity_t *p8est_connectivity_bcast (p8est_connectivity_t *
                                                conn_in, int root,
                                                sc_MPI_Comm comm);

/** Destroy a connectivity structure.  Also destroy all attributes.
 */
void                p8est_connectivity_destroy (p8est_connectivity_t *
                                                connectivity);

/** Allocate or free the attribute fields in a connectivity.
 * \param [in,out] conn         The conn->*_to_attr fields must either be NULL
 *                              or previously be allocated by this function.
 * \param [in] bytes_per_tree   If 0, tree_to_attr is freed (being NULL is ok).
 *                              If positive, requested space is allocated.
 */
void                p8est_connectivity_set_attr (p8est_connectivity_t * conn,
                                                 size_t bytes_per_tree);

/** Examine a connectivity structure.
 * \return  Returns true if structure is valid, false otherwise.
 */
int                 p8est_connectivity_is_valid (p8est_connectivity_t *
                                                 connectivity);

/** Check two connectivity structures for equality.
 * \return          Returns true if structures are equal, false otherwise.
 */
int                 p8est_connectivity_is_equal (p8est_connectivity_t * conn1,
                                                 p8est_connectivity_t *
                                                 conn2);

/** Write connectivity to a sink object.
 * \param [in] conn     The connectivity to be written.
 * \param [in,out] sink The connectivity is written into this sink.
 * \return              0 on success, nonzero on error.
 */
int                 p8est_connectivity_sink (p8est_connectivity_t * conn,
                                             sc_io_sink_t * sink);

/** Allocate memory and store the connectivity information there.
 * \param [in] conn     The connectivity structure to be exported to memory.
 * \param [in] code     Encoding and compression method for serialization.
 * \return              Newly created array that contains the information.
 */
sc_array_t         *p8est_connectivity_deflate (p8est_connectivity_t * conn,
                                                p8est_connectivity_encode_t
                                                code);

/** Save a connectivity structure to disk.
 * \param [in] filename         Name of the file to write.
 * \param [in] connectivity     Valid connectivity structure.
 * \return                      Returns 0 on success, nonzero on file error.
 */
int                 p8est_connectivity_save (const char *filename,
                                             p8est_connectivity_t *
                                             connectivity);

/** Read connectivity from a source object.
 * \param [in,out] source       The connectivity is read from this source.
 * \return              The newly created connectivity, or NULL on error.
 */
p8est_connectivity_t *p8est_connectivity_source (sc_io_source_t * source);

/** Create new connectivity from a memory buffer.
 * This function aborts on malloc errors.
 * \param [in] buffer   The connectivity is created from this memory buffer.
 * \return              The newly created connectivity, or NULL on format
 *                      error of the buffered connectivity data.
 */
p8est_connectivity_t *p8est_connectivity_inflate (sc_array_t * buffer);

/** Load a connectivity structure from disk.
 * \param [in] filename         Name of the file to read.
 * \param [out] bytes           Size in bytes of connectivity on disk or NULL.
 * \return              Returns valid connectivity, or NULL on file error.
 */
p8est_connectivity_t *p8est_connectivity_load (const char *filename,
                                               size_t *bytes);

/** Create a connectivity structure for the unit cube.
 */
p8est_connectivity_t *p8est_connectivity_new_unitcube (void);

/** Create a connectivity structure for an all-periodic unit cube.
 */
p8est_connectivity_t *p8est_connectivity_new_periodic (void);

/** Create a connectivity structure for a mostly periodic unit cube.
 * The left and right faces are identified, and bottom and top rotated.
 * Front and back are not identified.
 */
p8est_connectivity_t *p8est_connectivity_new_rotwrap (void);

/** Create a connectivity structure for a five-trees geometry with a
 * hole. The geometry is a 3D extrusion of the two drop example, and
 * covers [0, 3]*[0, 2]*[0, 3]. The additional dimension is Y.
 */
p8est_connectivity_t *p8est_connectivity_new_drop (void);

/** Create a connectivity structure that contains two cubes.
 */
p8est_connectivity_t *p8est_connectivity_new_twocubes (void);

/** Create a connectivity structure for two trees being rotated
 * w.r.t. each other in a user-defined way.
 * \param[in] l_face      index of left face
 * \param[in] r_face      index of right face
 * \param[in] orientation orientation of trees w.r.t. each other
 */
p8est_connectivity_t *p8est_connectivity_new_twotrees (int l_face,
                                                       int r_face,
                                                       int orientation);

/** Create a connectivity structure that contains two cubes
 * where the two far ends are identified periodically.
 */
p8est_connectivity_t *p8est_connectivity_new_twowrap (void);

/** Create a connectivity structure that contains a few cubes.
 * These are rotated against each other to stress the topology routines.
 */
p8est_connectivity_t *p8est_connectivity_new_rotcubes (void);

/** An m by n by p array with periodicity in x, y, and z if
 * periodic_a, periodic_b, and periodic_c are true, respectively.
 */
p8est_connectivity_t *p8est_connectivity_new_brick (int m, int n, int p,
                                                    int periodic_a,
                                                    int periodic_b,
                                                    int periodic_c);

/** Create a connectivity structure that builds a spherical shell.
 * It is made up of six connected parts [-1,1]x[-1,1]x[1,2].
 * This connectivity reuses vertices and relies on a geometry transformation.
 * It is thus not suitable for p8est_connectivity_complete.
 */
p8est_connectivity_t *p8est_connectivity_new_shell (void);

/** Create a connectivity structure that builds a solid sphere.
 * It is made up of two layers and a cube in the center.
 * This connectivity reuses vertices and relies on a geometry transformation.
 * It is thus not suitable for p8est_connectivity_complete.
 */
p8est_connectivity_t *p8est_connectivity_new_sphere (void);

/** Create a connectivity structure that builds a revolution torus.
 *
 * This connectivity reuses vertices and relies on a geometry transformation.
 * It is thus not suitable for p8est_connectivity_complete.
 *
 * This connectivity reuses ideas from disk2d connectivity.
 * More precisely the torus is divided into segments around
 * the revolution axis, each segments is made of 5 trees (à la disk2d).
 * The total number of trees if 5 times the number of segments.
 *
 * This connectivity is meant to be used with \ref p8est_geometry_new_torus
 *
 * \param[in] nSegments number of trees along the great circle
 */
p8est_connectivity_t *p8est_connectivity_new_torus (int nSegments);

/** Create connectivity structure from predefined catalogue.
 * \param [in]  name            Invokes connectivity_new_* function.
 *              brick235        brick (2, 3, 5, 0, 0, 0)
 *              periodic        periodic
 *              rotcubes        rotcubes
 *              rotwrap         rotwrap
 *              shell           shell
 *              sphere          sphere
 *              twocubes        twocubes
 *              twowrap         twowrap
 *              unit            unitcube
 * \return      An initialized connectivity if name is defined, NULL else.
 */
p8est_connectivity_t *p8est_connectivity_new_byname (const char *name);

/** Uniformly refine a connectivity.
 * This is useful if you would like to uniformly refine by something other
 * than a power of 2.
 *
 * \param [in] conn         A valid connectivity
 * \param [in] num_per_dim  The number of new trees in each direction.
 *                      Must use no more than \ref P8EST_OLD_QMAXLEVEL bits.
 *
 * \return a refined connectivity.
 */
p8est_connectivity_t *p8est_connectivity_refine (p8est_connectivity_t * conn,
                                                 int num_per_dim);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  iface       The number of the originating face.
 * \param [in]  nface       Encoded as nface = r * 6 + nf, where nf = 0..5 is
 *                          the neigbbor's connecting face number and r = 0..3
 *                          is the relative orientation to the neighbor's face.
 *                          This encoding matches p8est_connectivity_t.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0]..[2]    The coordinate axis sequence of the origin face,
 *                          the first two referring to the tangentials and the
 *                          third to the normal.  A permutation of (0, 1, 2).
 *              [3]..[5]    The coordinate axis sequence of the target face.
 *              [6]..[8]    Edge reversal flags for tangential axes (boolean);
 *                          face code in [0, 3] for the normal coordinate q:
 *                          0: q' = -q
 *                          1: q' = q + 1
 *                          2: q' = q - 1
 *                          3: q' = 2 - q
 */
void                p8est_expand_face_transform (int iface, int nface,
                                                 int ftransform[]);

/** Fill an array with the axis combination of a face neighbor transform.
 * \param [in]  connectivity Connectivity structure.
 * \param [in]  itree       The number of the originating tree.
 * \param [in]  iface       The number of the originating tree's face.
 * \param [out] ftransform  This array holds 9 integers.
 *              [0]..[2]    The coordinate axis sequence of the origin face.
 *              [3]..[5]    The coordinate axis sequence of the target face.
 *              [6]..[8]    Edge reversal flag for axes t1, t2; face code for n;
 *                          \see p8est_expand_face_transform.
 * \return                  The face neighbor tree if it exists, -1 otherwise.
 */
p4est_topidx_t      p8est_find_face_transform (p8est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iface, int ftransform[]);

/** Fills an array with information about edge neighbors.
 * \param [in] connectivity Connectivity structure.
 * \param [in] itree    The number of the originating tree.
 * \param [in] iedge    The number of the originating edge.
 * \param [in,out] ei   A p8est_edge_info_t structure with initialized array.
 */
void                p8est_find_edge_transform (p8est_connectivity_t *
                                               connectivity,
                                               p4est_topidx_t itree,
                                               int iedge,
                                               p8est_edge_info_t * ei);

/** Fills an array with information about corner neighbors.
 * \param [in] connectivity  Connectivity structure.
 * \param [in] itree    The number of the originating tree.
 * \param [in] icorner  The number of the originating corner.
 * \param [in,out] ci   A p8est_corner_info_t structure with initialized array.
 */
void                p8est_find_corner_transform (p8est_connectivity_t *
                                                 connectivity,
                                                 p4est_topidx_t itree,
                                                 int icorner,
                                                 p8est_corner_info_t * ci);

/** Internally connect a connectivity based on tree_to_vertex information.
 * Periodicity that is not inherent in the list of vertices will be lost.
 * \param [in,out] conn     The connectivity needs to have proper vertices
 *                          and tree_to_vertex fields.  The tree_to_tree
 *                          and tree_to_face fields must be allocated
 *                          and satisfy p8est_connectivity_is_valid (conn)
 *                          but will be overwritten.  The edge and corner
 *                          fields will be freed and allocated anew.
 */
void                p8est_connectivity_complete (p8est_connectivity_t * conn);

/** Removes corner and edge information of a connectivity
 *  such that enough information is left to run p8est_connectivity_complete successfully.
 *  The reduced connectivity still passes p8est_connectivity_is_valid.
 * \param [in,out] conn     The connectivity to be reduced.
 */
void                p8est_connectivity_reduce (p8est_connectivity_t * conn);

/** p8est_connectivity_permute
 * Given a permutation \a perm of the trees in a connectivity \a conn,
 * permute the trees of \a conn in place and update \a conn to match.
 * \param [in,out] conn                The connectivity whose trees are
 *                                     permuted.
 * \param [in] perm                    A permutation array, whose elements are
 *                                     size_t's.
 * \param [in] is_current_to_new       if true, the jth entry of perm is the
 *                                     new index for the entry whose current
 *                                     index is j, otherwise the jth entry of
 *                                     perm is the current index of the tree
 *                                     whose index will be j after the
 *                                     permutation.
 */
void                p8est_connectivity_permute (p8est_connectivity_t * conn,
                                                sc_array_t * perm,
                                                int is_current_to_new);

#ifdef P4EST_WITH_METIS

/** Reorder a connectivity using METIS.
 * This function takes a connectivity \a conn and a parameter \a k,
 * which will typically be the number of processes, and reorders the trees
 * such that if every processes is assigned (num_trees / k) trees, the
 * communication volume will be minimized.  This is intended for use with
 * connectivities that contain a large number of trees.  This should be done
 * BEFORE a p8est is created using the connectivity.  This is done in place:
 * any data structures that use indices to refer to trees before this
 * procedure will be invalid.  Note that this routine calls metis and not
 * parmetis because the connectivity is copied on every process.
 * A communicator is required because I'm not positive that metis is
 * deterministic. \a ctype determines when an edge exist between two trees in
 * the dual graph used by metis in the reordering.
 * \param [in]     comm       MPI communicator.
 * \param [in]     k          if k > 0, the number of pieces metis will use to
 *                            guide the reordering; if k = 0, the number of
 *                            pieces will be determined from the MPI
 *                            communicator.
 * \param [in,out] conn       connectivity that will be reordered.
 * \param [in]     ctype      determines when an edge exists in the dual graph
 *                            of the connectivity structure.
 */
void                p8est_connectivity_reorder (sc_MPI_Comm comm, int k,
                                                p8est_connectivity_t * conn,
                                                p8est_connect_type_t ctype);

/** Reorder a connectivity using METIS.
 * This is the same form of p8est_connectivity_reorder but it takes an initialized
 * sc array \a newid as extra argument.
 * In this way, the users can map old indices to new indices in the case it
 * is necessary (for instance to retrieve high-order nodes previously stored
 * in an array with old indices).
 * \param [in]     comm       MPI communicator.
 * \param [in]     k          if k > 0, the number of pieces metis will use to
 *                            guide the reordering; if k = 0, the number of
 *                            pieces will be determined from the MPI
 *                            communicator.
 * \param [in,out] conn       connectivity that will be reordered.
 * \param [in]     ctype      determines when an edge exists in the dual graph
 *                            of the connectivity structure.
 * \param [in,out] newid      array that maps old tree indices to new ones.
 *                            newid has to be an sc_array and it has to be
 *                            initialized (non-NULL) with element size
 *                            of size_t (using sc_array_new (sizeof (size_t))).
 *                            Input length arbitrary, output length modified.
 */
sc_array_t         *p8est_connectivity_reorder_newid (sc_MPI_Comm comm, int k,
                                                      p8est_connectivity_t *
                                                      conn,
                                                      p8est_connect_type_t
                                                      ctype,
                                                      sc_array_t * newid);

#endif /* P4EST_WITH_METIS */

/** p8est_connectivity_join_faces
 * This function takes an existing valid connectivity \a conn and modifies it
 * by joining two tree faces that are currently boundary faces.
 * \param [in,out] conn        connectivity that will be altered.
 * \param [in]     tree_left   tree that will be on the left side of the joined
 *                             faces.
 * \param [in]     tree_right  tree that will be on the right side of the
 *                             joined faces.
 * \param [in]     face_left   face of \a tree_left that will be joined.
 * \param [in]     face_right  face of \a tree_right that will be joined.
 * \param [in]     orientation the orientation of \a face_left and
 *                             \a face_right once joined (see the description
 *                             of p8est_connectivity_t to understand
 *                             orientation).
 */
void                p8est_connectivity_join_faces (p8est_connectivity_t *
                                                   conn,
                                                   p4est_topidx_t tree_left,
                                                   p4est_topidx_t tree_right,
                                                   int face_left,
                                                   int face_right,
                                                   int orientation);

/** p8est_connectivity_is_equivalent
 * This function compares two connectivities for equivalence: it returns
 * \a true if they are the same connectivity, or if they have the same
 * topology.  The definition of topological sameness is strict: there is no
 * attempt made to determine whether permutation and/or rotation of the trees
 * makes the connectivities equivalent.
 *
 * \param[in]      conn1    a valid connectivity
 * \param[out]     conn2    a valid connectivity
 */
int                 p8est_connectivity_is_equivalent (p8est_connectivity_t *
                                                      conn1,
                                                      p8est_connectivity_t *
                                                      conn2);

/** Return a pointer to a p8est_edge_transform_t array element. */
/*@unused@*/
static inline p8est_edge_transform_t *
p8est_edge_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_edge_transform_t));
  P4EST_ASSERT (it < array->elem_count);

  return (p8est_edge_transform_t *) (array->array +
                                     sizeof (p8est_edge_transform_t) * it);
}

/** Return a pointer to a p8est_corner_transform_t array element. */
/*@unused@*/
static inline p8est_corner_transform_t *
p8est_corner_array_index (sc_array_t * array, size_t it)
{
  P4EST_ASSERT (array->elem_size == sizeof (p8est_corner_transform_t));
  P4EST_ASSERT (it < array->elem_count);

  return
    (p8est_corner_transform_t *) (array->array +
                                  sizeof (p8est_corner_transform_t) * it);
}

/** Read an ABAQUS input file from a file stream.
 *
 * This utility function reads a basic ABAQUS file supporting element type with
 * the prefix C2D4, CPS4, and S4 in 2D and of type C3D8 reading them as
 * bilinear quadrilateral and trilinear hexahedral trees respectively.
 *
 * A basic 2D mesh is given below.  The \c *Node section gives the vertex
 * number and x, y, and z components for each vertex.  The \c *Element section
 * gives the 4 vertices in 2D (8 vertices in 3D) of each element in counter
 * clockwise order. So in 2D the nodes are given as:
 *
 *   4                     3
 *   +-------------------+
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   +-------------------+
 *   1                   2
 *
 * and in 3D they are given as:
 *
 * 8                     7
 *  +---------------------+
 *  |\                    |\
 *  | \                   | \
 *  |  \                  |  \
 *  |   \                 |   \
 *  |   5+---------------------+6
 *  |    |                |    |
 *  +----|----------------+    |
 *  4\   |               3 \   |
 *    \  |                  \  |
 *     \ |                   \ |
 *      \|                    \|
 *       +---------------------+
 *       1                     2
 *
 * \code
 * *Heading
 *  box.inp
 * *Node
 *     1,    5,   -5,    5
 *     2,    5,    5,    5
 *     3,    5,    0,    5
 *     4,   -5,    5,    5
 *     5,    0,    5,    5
 *     6,   -5,   -5,    5
 *     7,   -5,    0,    5
 *     8,    0,   -5,    5
 *     9,    0,    0,    5
 *    10,    5,    5,   -5
 *    11,    5,   -5,   -5
 *    12,    5,    0,   -5
 *    13,   -5,   -5,   -5
 *    14,    0,   -5,   -5
 *    15,   -5,    5,   -5
 *    16,   -5,    0,   -5
 *    17,    0,    5,   -5
 *    18,    0,    0,   -5
 *    19,   -5,   -5,    0
 *    20,    5,   -5,    0
 *    21,    0,   -5,    0
 *    22,   -5,    5,    0
 *    23,   -5,    0,    0
 *    24,    5,    5,    0
 *    25,    0,    5,    0
 *    26,    5,    0,    0
 *    27,    0,    0,    0
 * *Element, type=C3D8, ELSET=EB1
 *     1,       6,      19,      23,       7,       8,      21,      27,       9
 *     2,      19,      13,      16,      23,      21,      14,      18,      27
 *     3,       7,      23,      22,       4,       9,      27,      25,       5
 *     4,      23,      16,      15,      22,      27,      18,      17,      25
 *     5,       8,      21,      27,       9,       1,      20,      26,       3
 *     6,      21,      14,      18,      27,      20,      11,      12,      26
 *     7,       9,      27,      25,       5,       3,      26,      24,       2
 *     8,      27,      18,      17,      25,      26,      12,      10,      24
 * \endcode
 *
 * This code can be called two ways.  The first, when \c vertex==NULL and \c
 * tree_to_vertex==NULL, is used to count the number of trees and vertices in
 * the connectivity to be generated by the \c .inp mesh in the \a stream.  The
 * second, when \c vertices!=NULL and \c tree_to_vertex!=NULL, fill \c vertices
 * and \c tree_to_vertex.  In this case \c num_vertices and \c num_trees need
 * to be set to the maximum number of entries allocated in \c vertices and \c
 * tree_to_vertex.
 *
 * \param[in,out]  stream         file stream to read the connectivity from
 * \param[in,out]  num_vertices   the number of vertices in the connectivity
 * \param[in,out]  num_trees      the number of trees in the connectivity
 * \param[out]     vertices       the list of \c vertices of the connectivity
 * \param[out]     tree_to_vertex the \c tree_to_vertex map of the connectivity
 *
 * \returns 0 if successful and nonzero if not
 */

int                 p8est_connectivity_read_inp_stream (FILE * stream,
                                                        p4est_topidx_t *
                                                        num_vertices,
                                                        p4est_topidx_t *
                                                        num_trees,
                                                        double *vertices,
                                                        p4est_topidx_t *
                                                        tree_to_vertex);

/** Create a p4est connectivity from an ABAQUS input file.
 *
 * This utility function reads a basic ABAQUS file supporting element type with
 * the prefix C2D4, CPS4, and S4 in 2D and of type C3D8 reading them as
 * bilinear quadrilateral and trilinear hexahedral trees respectively.
 *
 * A basic 2D mesh is given below.  The \c *Node section gives the vertex
 * number and x, y, and z components for each vertex.  The \c *Element section
 * gives the 4 vertices in 2D (8 vertices in 3D) of each element in counter
 * clockwise order. So in 2D the nodes are given as:
 *
 *   4                     3
 *   +-------------------+
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   |                   |
 *   +-------------------+
 *   1                   2
 *
 * and in 3D they are given as:
 *
 * 8                     7
 *  +---------------------+
 *  |\                    |\
 *  | \                   | \
 *  |  \                  |  \
 *  |   \                 |   \
 *  |   5+---------------------+6
 *  |    |                |    |
 *  +----|----------------+    |
 *  4\   |               3 \   |
 *    \  |                  \  |
 *     \ |                   \ |
 *      \|                    \|
 *       +---------------------+
 *       1                     2
 *
 * \code
 * *Heading
 *  box.inp
 * *Node
 *     1,    5,   -5,    5
 *     2,    5,    5,    5
 *     3,    5,    0,    5
 *     4,   -5,    5,    5
 *     5,    0,    5,    5
 *     6,   -5,   -5,    5
 *     7,   -5,    0,    5
 *     8,    0,   -5,    5
 *     9,    0,    0,    5
 *    10,    5,    5,   -5
 *    11,    5,   -5,   -5
 *    12,    5,    0,   -5
 *    13,   -5,   -5,   -5
 *    14,    0,   -5,   -5
 *    15,   -5,    5,   -5
 *    16,   -5,    0,   -5
 *    17,    0,    5,   -5
 *    18,    0,    0,   -5
 *    19,   -5,   -5,    0
 *    20,    5,   -5,    0
 *    21,    0,   -5,    0
 *    22,   -5,    5,    0
 *    23,   -5,    0,    0
 *    24,    5,    5,    0
 *    25,    0,    5,    0
 *    26,    5,    0,    0
 *    27,    0,    0,    0
 * *Element, type=C3D8, ELSET=EB1
 *     1,       6,      19,      23,       7,       8,      21,      27,       9
 *     2,      19,      13,      16,      23,      21,      14,      18,      27
 *     3,       7,      23,      22,       4,       9,      27,      25,       5
 *     4,      23,      16,      15,      22,      27,      18,      17,      25
 *     5,       8,      21,      27,       9,       1,      20,      26,       3
 *     6,      21,      14,      18,      27,      20,      11,      12,      26
 *     7,       9,      27,      25,       5,       3,      26,      24,       2
 *     8,      27,      18,      17,      25,      26,      12,      10,      24
 * \endcode
 *
 * This function reads a mesh from \a filename and returns an associated p4est
 * connectivity.
 *
 * \param[in]  filename         file to read the connectivity from
 *
 * \returns an allocated connectivity associated with the mesh in \a filename
 */
p8est_connectivity_t *p8est_connectivity_read_inp (const char *filename);

SC_EXTERN_C_END;

#endif /* !P8EST_CONNECTIVITY_H */

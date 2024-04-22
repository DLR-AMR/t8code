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

/** \file p8est_mesh.h
 *
 * Forest topology in a conventional mesh format.
 *
 * A typical workflow starts with \ref p8est_mesh_params_init to initialize a
 * \ref p8est_mesh_params_t, followed by eventual user-dependent changes to the
 * parameters.
 *
 * Next a \ref p8est_mesh_t is created with \ref p8est_mesh_new_params.
 *
 * Now, the user can create a \ref p8est_mesh_face_neighbor_t with
 * \ref p8est_mesh_face_neighbor_init and loop over a quadrants face neighbors
 * by repeated calls to \ref p8est_mesh_face_neighbor_next.
 *
 * Once done, the mesh has to be destroyed with \ref p8est_mesh_destroy.
 *
 * \ingroup p8est
 */

#ifndef P8EST_MESH_H
#define P8EST_MESH_H

#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

/** This structure contains the different parameters of mesh creation.
 * A default instance can be initialzed by calling \ref p8est_mesh_params_init
 * and used for mesh creation by calling \ref p8est_mesh_new_params. */
typedef struct
{
  int                 compute_tree_index;     /**< Boolean to decide whether to
                                                   allocate and compute the
                                                   quad_to_tree list. */
  int                 compute_level_lists;    /**< Boolean to decide whether to
                                                   compute the level lists in
                                                   quad_level. */
  p8est_connect_type_t btype;                 /**< Flag indicating the
                                                   connection types (face, edge,
                                                   corner) stored in the mesh. */
  int                 edgehanging_corners;    /**< Boolean to decide whether to
                                                   add corner neighbors across
                                                   coarse edges. */
}
p8est_mesh_params_t;

/** This structure contains complete mesh information on a 2:1 balanced forest.
 * It stores the locally relevant neighborhood, that is, all locally owned
 * quadrants and one layer of adjacent ghost quadrants and their owners.
 *
 * For each local quadrant, its tree number is stored in quad_to_tree.
 * The quad_to_tree array is NULL by default and can be enabled using
 * \ref p8est_mesh_new_ext or \ref p8est_mesh_new_params.
 * For each ghost quadrant, its owner rank is stored in ghost_to_proc.
 * For each level, an array of local quadrant numbers is stored in quad_level.
 * The quad_level array is NULL by default and can be enabled using
 * \ref p8est_mesh_new_ext or \ref p8est_mesh_new_params.
 *
 * The quad_to_quad list stores one value for each local quadrant's face.
 * This value is in 0..local_num_quadrants-1 for local quadrants, or in
 * local_num_quadrants + (0..ghost_num_quadrants-1) for ghost quadrants.
 *
 * The quad_to_face list has equally many entries that are either:
 * 1. A value of v = 0..23 indicates one same-size neighbor.
 *    This value is decoded as v = r * 6 + nf, where nf = 0..5 is the
 *    neighbor's connecting face number and r = 0..3 is the relative
 *    orientation of the neighbor's face; see p8est_connectivity.h.
 * 2. A value of v = 24..119 indicates a double-size neighbor.
 *    This value is decoded as v = 24 + h * 24 + r * 6 + nf, where
 *    r and nf are as above and h = 0..3 is the number of the subface.
 *    h designates the subface of the large neighbor that the quadrant
 *    touches (this is the same as the large neighbor's face corner).
 * 3. A value of v = -24..-1 indicates four half-size neighbors.
 *    In this case the corresponding quad_to_quad index points into the
 *    quad_to_half array that stores four quadrant numbers per index,
 *    and the orientation of the smaller faces follows from 24 + v.
 *    The entries of quad_to_half encode between local and ghost quadrant
 *    in the same way as the quad_to_quad values described above.
 *    The small neighbors in quad_to_half are stored in the sequence
 *    of the face corners of this, i.e., the large quadrant.
 *
 * A quadrant on the boundary of the forest sees itself and its face number.
 *
 * The quad_to_edge list stores edge neighbors that are not face neighbors.
 * On the inside of a tree, there are one or two of those depending on size.
 * Between trees, there can be any number of same- or different-sized neighbors.
 * For same-tree same-size neighbors, we record their number in quad_to_edge
 * by the same convention as described for quad_to_quad above.  In this case,
 * the neighbor's matching edge number is always diagonally opposite,
 * that is, edge number ^ 3.
 *
 * For half- and double-size and all inter-tree edge neighbors, the
 * quad_to_edge value is in
 *    local_num_quadrants + local_num_ghosts + [0 .. local_num_edges - 1].
 * After subtracting the number of local and ghost quadrants,
 * it indexes into edge_offset, which encodes a group of edge neighbors.
 * Each member of a group may be one same/double-size quadrant or two half-size
 * quadrants; this is determined by the value of the edge_edge field as follows.
 * 1. A value of e = 0..23 indicates one same-size neighbor.
 *    This value is encoded as e = r * 12 + ne, where ne = 0..11 is the
 *    neighbor's connecting edge number and r = 0..1 indicates an edge flip.
 * 2. A value of e = 24..71 indicates a double-size neighbor.
 *    This value is decoded as e = 24 + h * 24 + r * 12 + ne, where
 *    r and ne are as above and h = 0..1 is the number of the subedge.
 *    h designates the subedge of the large neighbor that the quadrant
 *    touches (this is the same as the large neighbor's edge corner).
 * 3. A value of e = -24..-1 indicates two half-size neighbors.
 *    They are represented by two consecutive entries of the edge_quad and
 *    edge_edge arrays with identical values for edge_edge.
 *    The orientation of the smaller edges follows from 24 + e.
 *    The small neighbors in edge_quad are stored in the sequence
 *    of the edge corners of this, i.e., the large quadrant.
 *
 * Edges with no diagonal neighbor at all are assigned the value -3.  This
 * only happens on the domain boundary, which is necessarily a tree boundary.
 * Edge neighbors for face-hanging nodes are assigned the value -1.
 *
 * The quad_to_corner list stores corner neighbors that are not face or edge
 * neighbors.  On the inside of a tree, there is precisely one such neighbor
 * per corner.  In this case, its index is encoded as described above for
 * quad_to_quad.  The neighbor's matching corner number is always diagonally
 * opposite, that is, corner number ^ 7.
 *
 * On the inside of an inter-tree face, we have precisely one corner neighbor.
 * If a corner is across an inter-tree edge or corner, then the number of
 * corner neighbors may be any non-negative number.  In all three cases,
 * the quad_to_corner value is in
 *    local_num_quadrants + local_num_ghosts + [0 .. local_num_corners - 1].
 * After subtracting the number of local and ghost quadrants,
 * it indexes into corner_offset, which encodes a group of corner neighbors.
 * Each group contains the quadrant numbers encoded as usual for quad_to_quad
 * in corner_quad, and the corner number from the neighbor as corner_corner.
 *
 * Corners with no diagonal neighbor at all are assigned the value -3.  This
 * only happens on the domain boundary, which is necessarily a tree boundary.
 * If edgehanging_corners in params is 0, all corner-neighbors for face- and
 * edge-hanging nodes are assigned the value -1.
 * If it is 1, we check for corner neighbors across coarse edges and assign -1
 * for the remaining face- and edge-hanging nodes.
 *
 * The params struct describes the parameters the mesh was created with.
 * For full control over the parameters, use \ref p8est_mesh_new_params for
 * mesh creation.
 */
typedef struct
{
  p4est_locidx_t      local_num_quadrants; /**< number of process-local quadrants */
  p4est_locidx_t      ghost_num_quadrants; /**< number of ghost-layer quadrants */

  p4est_topidx_t     *quad_to_tree;     /**< tree index for each local quad.
                                             Is NULL if compute_tree_index in
                                             params is 0. */
  int                *ghost_to_proc;    /**< processor for each ghost quad */

  p4est_locidx_t     *quad_to_quad;     /**< one index for each of the 6 faces */
  int8_t             *quad_to_face;     /**< encodes orientation/2:1 status */
  sc_array_t         *quad_to_half;     /**< stores half-size neighbors */
  sc_array_t         *quad_level;       /**< Stores lists of per-level quads.
                                             The array has entries indexed by
                                             0..P4EST_QMAXLEVEL inclusive that
                                             are arrays of local quadrant ids.
                                             Is NULL if compute_level_lists in
                                             params is 0. */

  /* These members are NULL if btype in params is < P8EST_CONNECT_EDGE and can
   * be requested in \ref p8est_mesh_new. */
  p4est_locidx_t      local_num_edges;  /**< unsame-size and tree-boundary edges */
  p4est_locidx_t     *quad_to_edge;     /**< 12 indices for each local quad */
  sc_array_t         *edge_offset;      /**< local_num_edges + 1 entries */
  sc_array_t         *edge_quad;        /**< edge_offset indexes into this */
  sc_array_t         *edge_edge;        /**< and this one too (type int8_t) */

  /* These members are NULL if btype in params is < P8EST_CONNECT_CORNER and can
   * be requested in \ref p8est_mesh_new. */
  p4est_locidx_t      local_num_corners;        /**< tree-boundary corners */
  p4est_locidx_t     *quad_to_corner;   /**< 8 indices for each local quad */
  sc_array_t         *corner_offset;    /**< local_num_corners + 1 entries */
  sc_array_t         *corner_quad;      /**< corner_offset indexes into this */
  sc_array_t         *corner_corner;    /**< and this one too (type int8_t) */

  p8est_mesh_params_t params;           /**< parameters the mesh was created
                                             with, e.g. by passing them to
                                             \ref p8est_mesh_new_ext or
                                             \ref p8est_mesh_new_params */
}
p8est_mesh_t;

/** This structure can be used as the status of a face neighbor iterator.
  * It always contains the face and subface of the neighbor to be processed.
  */
typedef struct
{
  /* forest information */
  p8est_t            *p4est;            /**< the forest */
  p8est_ghost_t      *ghost;            /**< the ghost layer of the forest */
  p8est_mesh_t       *mesh;             /**< a mesh derived from the forest */

  /* quadrant information */
  p4est_topidx_t      which_tree;       /**< the current tree index */
  p4est_locidx_t      quadrant_id;      /**< tree-local quadrant index */
  p4est_locidx_t      quadrant_code;    /**< 6 * (quadrant_id + tree_offset) */

  /* neighbor information */
  int                 face;     /**< Face number in 0..5. */
  int                 subface;  /**< Hanging neighbor number in 0..3. */

  /* internal information */
  p4est_locidx_t      current_qtq;      /**< track index of current neighboring
                                             quadrant */
}
p8est_mesh_face_neighbor_t;

/** Calculate the memory usage of the mesh structure.
 * \param [in] mesh     Mesh structure.
 * \return              Memory used in bytes.
 */
size_t              p8est_mesh_memory_used (p8est_mesh_t * mesh);

/** Initialize a default p8est_mesh_params_t structure.
 * The parameters are set to create the most basic mesh structure, without
 * tree index and level lists and considering only face connections. */
void                p8est_mesh_params_init (p8est_mesh_params_t * params);

/** Create a p8est_mesh structure.
 * This function does not populate the quad_to_tree and quad_level fields and
 * ignores corner neighbors across edge-hanging corners.
 * To populate them, use \ref p8est_mesh_new_params.
 * \param [in] p8est    A forest that is fully 2:1 balanced.
 * \param [in] ghost    The ghost layer created from the provided p4est.
 * \param [in] btype    Determines the highest codimension of neighbors.
 * \return              A fully allocated mesh structure.
 */
p8est_mesh_t       *p8est_mesh_new (p8est_t * p8est, p8est_ghost_t * ghost,
                                    p8est_connect_type_t btype);

/** Create a new mesh.
 * \param [in] p8est    A forest that is fully 2:1 balanced.
 * \param [in] ghost    The ghost layer created from the provided p4est.
 * \param [in] params   The mesh creation parameters. If NULL, the function
 *                      defaults to the parameters of \ref p8est_mesh_params_init.
 * \return              A fully allocated mesh structure.
 */
p8est_mesh_t       *p8est_mesh_new_params (p8est_t * p8est,
                                           p8est_ghost_t * ghost,
                                           p8est_mesh_params_t * params);

/** Destroy a p8est_mesh structure.
 * \param [in] mesh     Mesh structure previously created by p8est_mesh_new.
 */
void                p8est_mesh_destroy (p8est_mesh_t * mesh);

/** Access a process-local quadrant inside a forest.
 * Needs a mesh with populated quad_to_tree array.
 * This is a special case of \ref p8est_mesh_quadrant_cumulative.
 *
 * \param [in] p4est  The forest.
 * \param [in] mesh   The mesh.
 * \param [in] qid    Process-local id of the quadrant (cumulative over trees).
 * \return            A pointer to the requested quadrant.
 */
p8est_quadrant_t   *p8est_mesh_get_quadrant (p8est_t * p4est,
                                             p8est_mesh_t * mesh,
                                             p4est_locidx_t qid);

/*** OUTDATED FUNCTION ***/
/** Lookup neighboring quads of quadrant in a specific direction
 * \param [in]  p4est              Forest to be worked with.
 * \param [in]  ghost              Ghost quadrants.
 * \param [in]  mesh               Mesh structure.
 * \param [in]  curr_quad_id       Process-local ID of current quad.
 * \param [in]  direction          Direction in which to look for adjacent
 *                                 quadrants is encoded as follows:
 *                                   0 ..  5 neighbor(-s) across f_i,
 *                                   6 .. 17 neighbor(-s) across e_{i-6}
 *                                  18 .. 25 neighbor(-s) across c_{i-18}
 * \param [out] neighboring_quads  Array containing neighboring quad(-s)
 *                                 Needs to be empty, contains
 *                                 p4est_quadrant_t*. May be NULL, then \b
 *                                 neighboring_qids must not be NULL.
 * \param [out] neighboring_encs   Array containing encodings for neighboring
 *                                 quads as described below
 *                                 Needs to be empty, contains int.
 * CAUTION: Note, that the encodings differ from the encodings saved in the
 *          mesh.
 *          Positive values are for local quadrants, negative values indicate
 *          ghost quadrants.
 *          Faces:     1 ..  24 => same size neighbor
 *                                 (r * 6 + nf) + 1; nf = 0 .. 5 face index;
 *                                 r = 0 .. 3 relative orientation
 *                    25 .. 120 => double size neighbor
 *                                 25 + h * 24 + r * 6 + nf; h = 0 .. 3 number
 *                                 of the subface; r, nf as above
 *                   121 .. 144 => half size neighbors
 *                                 121 + r * 6 + nf; r, nf as above
 *          Edges:     1 ..  24 => same size neighbor
 *                                 r * 12 + ne + 1; ne = 0 .. 11 edge index;
 *                                 r = 0 .. 1 relative orientation
 *                    25 ..  72 => double size neighbor
 *                                 25 + h * 24 + r * 12 + ne; h = 0 .. 1 number
 *                                 of the subedge; r, ne as above
 *                    73 ..  96 => half size neighbors
 *                                 73 + r * 12 + ne; r, ne as above
 *          Corners:   1 ..   8 => nc + 1; nc = 0 .. 7 corner index
 * \param [out] neighboring_qids   Array containing quadrant ids for neighboring
 *                                 quadrants. May be NULL, then no neighboring
 *                                 qids are collected.
 *                                 If non-NULL the array needs to be empty and
 *                                 will contain int.
 */
p4est_locidx_t      p8est_mesh_get_neighbors (p8est_t * p4est,
                                              p8est_ghost_t * ghost,
                                              p8est_mesh_t * mesh,
                                              p4est_locidx_t curr_quad_id,
                                              p4est_locidx_t direction,
                                              sc_array_t * neighboring_quads,
                                              sc_array_t * neighboring_encs,
                                              sc_array_t * neighboring_qids);

/** Find a quadrant based on its cumulative number in the local forest.
 * If the quad_to_tree field of the mesh structure exists, this is O(1).
 * Otherwise, we perform a binary search over the processor-local trees.
 *
 * \param [in]  p8est           Forest to be worked with.
 * \param [in]  mesh            A mesh derived from the forest.
 * \param [in]  cumulative_id   Cumulative index over all trees of quadrant.
 *                              Must refer to a local (non-ghost) quadrant.
 * \param [in,out] which_tree   If not NULL, the input value can be -1
 *                              or an initial guess for the quadrant's tree
 *                              and output is the tree of returned quadrant.
 * \param [out] quadrant_id     If not NULL, the number of quadrant in tree.
 * \return                      The identified quadrant.
 */
p8est_quadrant_t   *p8est_mesh_quadrant_cumulative (p8est_t * p8est,
                                                    p8est_mesh_t * mesh,
                                                    p4est_locidx_t
                                                    cumulative_id,
                                                    p4est_topidx_t *
                                                    which_tree,
                                                    p4est_locidx_t *
                                                    quadrant_id);

/** Initialize a mesh neighbor iterator by quadrant index.
 * \param [out] mfn         A p8est_mesh_face_neighbor_t to be initialized.
 * \param [in]  p8est       Forest to be worked with.
 * \param [in]  ghost       Ghost layer of the forest.
 * \param [in]  mesh        A mesh derived from the forest.
 * \param [in]  which_tree  Tree of quadrant whose neighbors are looped over.
 * \param [in]  quadrant_id Index relative to which_tree of quadrant.
 */
void                p8est_mesh_face_neighbor_init2 (p8est_mesh_face_neighbor_t
                                                    * mfn, p8est_t * p8est,
                                                    p8est_ghost_t * ghost,
                                                    p8est_mesh_t * mesh,
                                                    p4est_topidx_t which_tree,
                                                    p4est_locidx_t
                                                    quadrant_id);

/** Initialize a mesh neighbor iterator by quadrant pointer.
 * \param [out] mfn         A p8est_mesh_face_neighbor_t to be initialized.
 * \param [in]  p8est       Forest to be worked with.
 * \param [in]  ghost       Ghost layer of the forest.
 * \param [in]  mesh        A mesh derived from the forest.
 * \param [in]  which_tree  Tree of quadrant whose neighbors are looped over.
 * \param [in]  quadrant    Pointer to quadrant contained in which_tree.
 */
void                p8est_mesh_face_neighbor_init (p8est_mesh_face_neighbor_t
                                                   * mfn, p8est_t * p8est,
                                                   p8est_ghost_t * ghost,
                                                   p8est_mesh_t * mesh,
                                                   p4est_topidx_t which_tree,
                                                   p8est_quadrant_t *
                                                   quadrant);

/** Move the iterator forward to loop around neighbors of the quadrant.
 * \param [in,out] mfn      Internal status of the iterator.
 * \param [out]    ntree    If not NULL, the tree number of the neighbor.
 * \param [out]    nquad    If not NULL, the quadrant number within tree.
 *                          For ghosts instead the number in ghost layer.
 * \param [out]    nface    If not NULL, neighbor's face encoding as in
                            quad_to_face array of p8est_mesh_t.
 * \param [out]    nrank    If not NULL, the owner process of the neighbor.
 * \return                  Either a real quadrant or one from the ghost layer.
 *                          Returns NULL when the iterator is done.
 */
p8est_quadrant_t   *p8est_mesh_face_neighbor_next (p8est_mesh_face_neighbor_t
                                                   * mfn,
                                                   p4est_topidx_t * ntree,
                                                   p4est_locidx_t * nquad,
                                                   int *nface, int *nrank);

/** Get the user data for the current face neighbor.
 * \param [in]     mfn           Internal status of the iterator.
 * \param [in]     ghost_data    Data for the ghost quadrants that has been
 *                               synchronized with p4est_ghost_exchange_data.
 * \return                       A pointer to the user data for the current
 *                               neighbor.
 */
void               *p8est_mesh_face_neighbor_data (p8est_mesh_face_neighbor_t
                                                   * mfn, void *ghost_data);

SC_EXTERN_C_END;

#endif /* !P8EST_MESH_H */

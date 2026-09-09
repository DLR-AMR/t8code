/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

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

/** \file t8_cmesh_helpers.hxx
 *
 * Shared coordinate-matching primitives used to find tree vertices that
 * coincide geometrically (within a tolerance derived from the mesh's
 * bounding box), via a spatial hash ("bin") of quantized coordinates, and
 * \ref t8_cmesh_join_by_vertices, which uses them to derive face
 * connectivity of an uncommitted cmesh from its trees' corner coordinates
 * or global vertex ids.
 */

#ifndef T8_CMESH_HELPERS_HXX
#define T8_CMESH_HELPERS_HXX

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_types/t8_vec.hxx>

#include <array>
#include <vector>

/** Holds the domain bounds, tolerance and bin scale shared by all
 * coordinate-hash and coordinate-match operations for one mesh.
 */
struct t8_cmesh_coord_match_context
{
  double min_coord;        /**< Minimum coordinate value over all given tree vertices. */
  double max_coord;        /**< Maximum coordinate value over all given tree vertices. */
  double tolerance;        /**< Scaled tolerance to decide if two doubles are equal. */
  double inverse_bin_size; /**< Number of bins divided by the domain's coordinate range. */
};

/** Scan all given tree vertices for their coordinate bounds and derive a
 * scaled tolerance and bin scale from them.
 * \param [in] ntrees    Number of coarse mesh elements resp. trees.
 * \param [in] eclasses  Element classes of the trees, of length \a ntrees.
 * \param [in] vertices  Per-tree vertex coordinates, one corner array per tree.
 * \return The resulting coordinate match context.
 */
t8_cmesh_coord_match_context
t8_cmesh_compute_coord_match_context (const t8_gloidx_t ntrees, const t8_eclass_t *eclasses,
                                      const std::vector<std::vector<t8_3D_vec>> &vertices);

/** Compute a spatial-hash bin key for a single tree vertex's coordinates.
 * Vertices at (nearly) the same location hash to the same key, so summing
 * this over a face's vertices reproduces a face-level hash with the same
 * property.
 * \param [in] ctx     A context as computed by \ref t8_cmesh_compute_coord_match_context.
 * \param [in] vertex  The vertex's coordinates.
 * \return The hash key.
 */
unsigned long
t8_cmesh_hash_coord (const t8_cmesh_coord_match_context &ctx, const t8_3D_vec &vertex);

/** Check whether two tree vertices coincide within \a ctx's tolerance.
 * \param [in] ctx       A context as computed by \ref t8_cmesh_compute_coord_match_context.
 * \param [in] vertex_a  The first vertex's coordinates.
 * \param [in] vertex_b  The second vertex's coordinates.
 * \return True if \a vertex_a and \a vertex_b agree within \a ctx's tolerance in every component.
 */
bool
t8_cmesh_coords_match (const t8_cmesh_coord_match_context &ctx, const t8_3D_vec &vertex_a, const t8_3D_vec &vertex_b);

/** Set the face connectivity of \a cmesh by matching tree vertices, using
 * \a cmesh's stash. Two faces are joined if all their vertices match.
 * Vertices are compared by their already-assigned global vertex id if
 * \a cmesh's stash has any (see \ref t8_cmesh_set_global_vertices_of_tree,
 * \ref t8_cmesh_vertex_conn_set_vertices_by_coordinates); otherwise ids are
 * derived from coordinates internally, purely to speed up the matching, and
 * are discarded afterwards. Either way, this routine never enables or
 * changes \a cmesh's vertex connectivity itself: callers who also want
 * global vertex ids on the committed cmesh have to request that
 * independently (\ref t8_cmesh_enable_vertex_conn /
 * \ref t8_cmesh_vertex_conn_set_vertices_by_coordinates).
 * Faces that already have an explicit join (\ref t8_cmesh_set_join) are
 * left untouched.
 *
 * \param [in,out] cmesh  An initialized, uncommitted cmesh whose stash still
 *                        holds tree classes and, depending on the matching
 *                        strategy used, vertex coordinates or global vertex ids.
 *
 * \note Must be called before the stash is converted into the committed
 *       tree structure (i.e. before \a cmesh's stash is consumed/destroyed).
 * \note This routine does not detect periodic boundaries.
 * \note This routine might be too expensive for very large meshes if it
 *       has to derive vertex ids from coordinates. In this case, consider
 *       to use a fully featured mesh generator.
 */
void
t8_cmesh_join_by_vertices (const t8_cmesh_t cmesh);

#endif /* !T8_CMESH_HELPERS_HXX */

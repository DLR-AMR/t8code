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

#include <t8_forest/t8_forest_touching_ranks.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <algorithm>
#include <iterator>
#include <set>
#include <utility>

/** A neighbor of a cover element, routed to the ranks that own a part of it. */
struct t8_forest_touching_query
{
  t8_forest_position first;        /**< The position of the first descendant of the neighbor. */
  t8_forest_position last;         /**< The position of the last descendant of the neighbor. */
  t8_linearidx_t id;               /**< The linear id of the neighbor at its level. */
  int level;                       /**< The level of the neighbor. */
  int bdy_dim;                     /**< The dimension of the boundary of the neighbor shared with the cover element. */
  int bdy_id;                      /**< The id of this boundary of the neighbor. */
  int source;                      /**< The rank of the cover element. */
  t8_forest_position source_first; /**< The first position of the source rank. */
  t8_forest_position source_next;  /**< The first position after the source rank. */

  /** Order by all members, to remove duplicates. */
  friend auto
  operator<=> (const t8_forest_touching_query &, const t8_forest_touching_query &)
    = default;
};

/** The maximal elements between the first position of this rank and the first position of the next rank with
 * elements, in SFC order, each with its local tree. The caller destroys the elements.
 * \param [in] forest  A committed forest with splitters.
 * \param [out] cover  The cover elements are appended.
 */
static void
t8_forest_partition_cover (const t8_forest_t forest, std::vector<std::pair<t8_locidx_t, t8_element_t *>> &cover)
{
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  if (!splitters.has_elements) {
    return;
  }
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const int maxlevel = forest->maxlevel;
  t8_forest_position position = splitters.first_position;
  while (position < splitters.next_position) {
    const t8_locidx_t ltreeid = position.gtreeid - forest->first_local_tree;
    T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
    t8_element_t *element, *parent, *desc;
    scheme->element_new (eclass, 1, &element);
    scheme->element_new (eclass, 1, &parent);
    scheme->element_new (eclass, 1, &desc);
    /* Coarsen the descendant at the position as long as the parent starts at the position and lies in the range */
    scheme->element_set_linear_id (eclass, element, maxlevel, position.id);
    while (scheme->element_get_level (eclass, element) > 0) {
      scheme->element_get_parent (eclass, element, parent);
      scheme->element_get_first_descendant (eclass, parent, desc, maxlevel);
      if (scheme->element_get_linear_id (eclass, desc, maxlevel) != position.id) {
        break;
      }
      scheme->element_get_last_descendant (eclass, parent, desc, maxlevel);
      const t8_forest_position parent_last
        = { position.gtreeid, scheme->element_get_linear_id (eclass, desc, maxlevel) };
      if (!(parent_last < splitters.next_position)) {
        break;
      }
      scheme->element_copy (eclass, parent, element);
    }
    cover.emplace_back (ltreeid, element);
    /* Continue after the last descendant of the element, in the next tree if it is the last one of its tree */
    scheme->element_get_last_descendant (eclass, element, desc, maxlevel);
    const t8_linearidx_t last_id = scheme->element_get_linear_id (eclass, desc, maxlevel);
    scheme->element_set_linear_id (eclass, parent, 0, 0);
    scheme->element_get_last_descendant (eclass, parent, desc, maxlevel);
    if (last_id == scheme->element_get_linear_id (eclass, desc, maxlevel)) {
      position = { position.gtreeid + 1, 0 };
    }
    else {
      position = { position.gtreeid, last_id + 1 };
    }
    scheme->element_destroy (eclass, 1, &parent);
    scheme->element_destroy (eclass, 1, &desc);
  }
}

/** Whether the part of this rank touches a boundary of an element with dimension at least min_dim. If the element
 * lies partly in the part, the children touching the boundary are checked (the children at the face if only faces
 * are considered).
 * \param [in] forest   A committed forest with splitters.
 * \param [in] gtreeid  The global tree of \a element, a local tree.
 * \param [in] eclass   The eclass of the tree.
 * \param [in] element  An element of the tree.
 * \param [in] bdy_dim  The dimension of the boundary of \a element.
 * \param [in] bdy_id   The id of the boundary.
 * \param [in] min_dim  The minimal dimension of the contact, faces only if it is the dimension of the tree minus one.
 * \return True if the part of this rank touches the boundary with dimension at least \a min_dim.
 */
static bool
t8_forest_part_touches_boundary (const t8_forest_t forest, const t8_gloidx_t gtreeid, const t8_eclass_t eclass,
                                 const t8_element_t *element, const int bdy_dim, const int bdy_id, const int min_dim)
{
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  t8_forest_position first, last;
  t8_forest_element_positions (forest, gtreeid, eclass, element, &first, &last);
  if (last < splitters.first_position || !(first < splitters.next_position)) {
    return false;
  }
  if (splitters.first_position <= first && last < splitters.next_position) {
    return true;
  }
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  T8_ASSERT (scheme->element_get_level (eclass, element) < forest->maxlevel);
  bool touches = false;
  if (min_dim == t8_eclass_to_dimension[eclass] - 1) {
    const int num_children = scheme->element_get_num_face_children (eclass, element, bdy_id);
    std::vector<t8_element_t *> children (num_children);
    scheme->element_new (eclass, num_children, children.data ());
    scheme->element_get_children_at_face (eclass, element, bdy_id, children.data (), num_children, NULL);
    for (int ichild = 0; ichild < num_children && !touches; ichild++) {
      const int child_face = scheme->element_face_get_child_face (eclass, element, bdy_id, ichild);
      touches
        = t8_forest_part_touches_boundary (forest, gtreeid, eclass, children[ichild], bdy_dim, child_face, min_dim);
    }
    scheme->element_destroy (eclass, num_children, children.data ());
  }
  else {
    const int num_children = scheme->element_get_num_boundary_touching_children (eclass, element, bdy_dim, bdy_id);
    std::vector<t8_element_t *> children (num_children);
    std::vector<int> child_bdy_dims (num_children), child_bdy_ids (num_children);
    scheme->element_new (eclass, num_children, children.data ());
    scheme->element_get_children_touching_boundary (eclass, element, bdy_dim, bdy_id, children.data (), num_children,
                                                    child_bdy_dims.data (), child_bdy_ids.data ());
    for (int ichild = 0; ichild < num_children && !touches; ichild++) {
      if (child_bdy_dims[ichild] >= min_dim) {
        touches = t8_forest_part_touches_boundary (forest, gtreeid, eclass, children[ichild], child_bdy_dims[ichild],
                                                   child_bdy_ids[ichild], min_dim);
      }
    }
    scheme->element_destroy (eclass, num_children, children.data ());
  }
  return touches;
}

void
t8_forest_touching_ranks_compute (t8_forest_t forest, const int min_dim)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const int dim = forest->dimension;
  const bool faces_only = min_dim == dim - 1;
  T8_ASSERT (0 <= min_dim && min_dim <= SC_MAX (dim - 1, 0));
  if (!faces_only) {
    /* The neighbors across tree boundaries need the vertex and edge connectivity of the trees */
    t8_cmesh_build_boundary_connectivity (t8_forest_get_cmesh (forest), forest->mpicomm);
  }

  /* The neighbors of the cover elements that do not lie in the own range */
  std::vector<t8_forest_touching_query> queries;
  const auto add_query = [&] (const t8_gloidx_t gtreeid, const t8_eclass_t eclass, const t8_element_t *neigh,
                              const int bdy_dim, const int bdy_id) {
    t8_forest_touching_query query;
    t8_forest_element_positions (forest, gtreeid, eclass, neigh, &query.first, &query.last);
    if (splitters.first_position <= query.first && query.last < splitters.next_position) {
      return;
    }
    query.level = scheme->element_get_level (eclass, neigh);
    query.id = scheme->element_get_linear_id (eclass, neigh, query.level);
    query.bdy_dim = bdy_dim;
    query.bdy_id = bdy_id;
    query.source = forest->mpirank;
    query.source_first = splitters.first_position;
    query.source_next = splitters.next_position;
    queries.push_back (query);
  };
  std::vector<std::pair<t8_locidx_t, t8_element_t *>> cover;
  t8_forest_partition_cover (forest, cover);
  sc_array_t neighbors;
  sc_array_init (&neighbors, sizeof (t8_forest_boundary_neighbor_t));
  for (auto &[ltreeid, element] : cover) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
    if (faces_only) {
      const int num_faces = scheme->element_get_num_faces (eclass, element);
      for (int face = 0; face < num_faces; face++) {
        const t8_eclass_t neigh_eclass = t8_forest_element_neighbor_eclass (forest, ltreeid, element, face);
        t8_element_t *neigh;
        scheme->element_new (neigh_eclass, 1, &neigh);
        int dual_face;
        const t8_gloidx_t neigh_tree
          = t8_forest_element_face_neighbor (forest, ltreeid, element, neigh, neigh_eclass, face, &dual_face);
        if (neigh_tree >= 0) {
          add_query (neigh_tree, neigh_eclass, neigh, dim - 1, dual_face);
        }
        scheme->element_destroy (neigh_eclass, 1, &neigh);
      }
    }
    else {
      for (int bdy_dim = dim - 1; bdy_dim >= min_dim; bdy_dim--) {
        const int num_boundaries = scheme->element_get_num_boundaries (eclass, element, bdy_dim);
        for (int bdy_id = 0; bdy_id < num_boundaries; bdy_id++) {
          t8_forest_element_boundary_neighbors (forest, ltreeid, element, bdy_dim, bdy_id, &neighbors);
          for (size_t ineigh = 0; ineigh < neighbors.elem_count; ineigh++) {
            const t8_forest_boundary_neighbor_t *neighbor
              = (const t8_forest_boundary_neighbor_t *) sc_array_index (&neighbors, ineigh);
            add_query (neighbor->gtreeid, neighbor->eclass, neighbor->neigh, bdy_dim, neighbor->bdy_id);
          }
          t8_forest_boundary_neighbors_clear (forest, &neighbors);
        }
      }
    }
    scheme->element_destroy (eclass, 1, &element);
  }
  sc_array_reset (&neighbors);
  /* Cover elements of the same level may share a neighbor at the same boundary of it */
  std::sort (queries.begin (), queries.end ());
  queries.erase (std::unique (queries.begin (), queries.end ()), queries.end ());

  /* Route each neighbor to the ranks that own a part of it */
  t8_forest_route (forest, queries, [] (const t8_forest_touching_query &query, const t8_forest_bisection &level) {
    return t8_forest_route_halves_positions (level, query.first, query.last, false);
  });

  /* The senders whose cover element shares a boundary with the part of this rank are the touching ranks */
  t8_forest_touching_ranks_destroy (forest);
  t8_forest_touching_ranks *touching = forest->touching_ranks = new t8_forest_touching_ranks;
  touching->min_dim = min_dim;
  if (splitters.has_elements) {
    touching->ranges.push_back ({ forest->mpirank, splitters.first_position, splitters.next_position });
  }
  std::set<int> touching_sources;
  for (const t8_forest_touching_query &query : queries) {
    if (query.source == forest->mpirank || touching_sources.count (query.source) > 0) {
      continue;
    }
    const t8_locidx_t ltreeid = query.first.gtreeid - forest->first_local_tree;
    T8_ASSERT (0 <= ltreeid && ltreeid < t8_forest_get_num_local_trees (forest));
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
    t8_element_t *neigh;
    scheme->element_new (eclass, 1, &neigh);
    scheme->element_set_linear_id (eclass, neigh, query.level, query.id);
    if (t8_forest_part_touches_boundary (forest, query.first.gtreeid, eclass, neigh, query.bdy_dim, query.bdy_id,
                                         min_dim)) {
      touching->ranges.push_back ({ query.source, query.source_first, query.source_next });
      touching_sources.insert (query.source);
    }
    scheme->element_destroy (eclass, 1, &neigh);
  }
  std::sort (touching->ranges.begin (), touching->ranges.end (),
             [] (const t8_forest_rank_range &a, const t8_forest_rank_range &b) { return a.rank < b.rank; });
}

void
t8_forest_touching_ranks_require (t8_forest_t forest, const int min_dim)
{
  if (forest->touching_ranks == NULL || forest->touching_ranks->min_dim > min_dim) {
    t8_forest_touching_ranks_compute (forest, min_dim);
  }
}

void
t8_forest_touching_ranks_copy (t8_forest_t forest, const t8_forest_t from)
{
  t8_forest_touching_ranks_destroy (forest);
  if (from->touching_ranks != NULL) {
    forest->touching_ranks = new t8_forest_touching_ranks (*from->touching_ranks);
  }
}

void
t8_forest_touching_ranks_destroy (t8_forest_t forest)
{
  delete forest->touching_ranks;
  forest->touching_ranks = NULL;
}

int
t8_forest_touching_ranks_match_splitters (const t8_forest_t forest)
{
  T8_ASSERT (forest->touching_ranks != NULL);
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  for (const t8_forest_rank_range &range : forest->touching_ranks->ranges) {
    if (range.rank == forest->mpirank) {
      return splitters.has_elements && range.first == splitters.first_position && range.next == splitters.next_position;
    }
  }
  return !splitters.has_elements;
}

void
t8_forest_touching_ranks_owner_bounds (const t8_forest_t forest, const t8_forest_position &position, int *lower,
                                       int *upper)
{
  T8_ASSERT (forest->touching_ranks != NULL);
  const std::vector<t8_forest_rank_range> &ranges = forest->touching_ranks->ranges;
  /* The first range that starts after the position */
  const auto after = std::upper_bound (
    ranges.begin (), ranges.end (), position,
    [] (const t8_forest_position &pos, const t8_forest_rank_range &range) { return pos < range.first; });
  if (after != ranges.begin () && position < std::prev (after)->next) {
    *lower = *upper = std::prev (after)->rank;
    return;
  }
  /* The owner lies between the touching ranks before and after the position */
  *lower = after == ranges.begin () ? 0 : std::prev (after)->rank + 1;
  *upper = after == ranges.end () ? forest->mpisize - 1 : after->rank - 1;
  T8_ASSERT (*lower <= *upper);
}

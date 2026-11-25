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

#include <t8_forest/t8_forest_ghost_overlap.hxx>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_cmesh/t8_cmesh.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <algorithm>
#include <limits>
#include <set>
#include <tuple>

/** An element of the overlap region of a rank, routed to the ranks that own a part of it. */
struct t8_forest_overlap_query
{
  t8_forest_position first; /**< The position of the first descendant of the element. */
  t8_forest_position last;  /**< The position of the last descendant of the element. */
  t8_linearidx_t id;        /**< The linear id of the element at its level. */
  int level;                /**< The level of the element. */
  int source;               /**< The rank whose region contains the element. */
  t8_forest_position source_first; /**< The first position of the source rank, to send the notice to. */

  /** Order by all members, to remove duplicates. */
  friend auto
  operator<=> (const t8_forest_overlap_query &, const t8_forest_overlap_query &)
    = default;
};

/** A notice that a rank sends overlap ghosts to the owner of a position. */
struct t8_forest_overlap_notice
{
  t8_forest_position target; /**< The first position of the receiving rank. */
  int sender;                /**< The rank that sends overlap ghosts to it. */

  /** Order by all members, to remove duplicates. */
  friend auto
  operator<=> (const t8_forest_overlap_notice &, const t8_forest_overlap_notice &)
    = default;
};

/** An element of a region, identified by its tree, level and linear id. */
struct t8_forest_region_element
{
  t8_gloidx_t gtreeid; /**< The global tree. */
  int level;           /**< The level. */
  t8_linearidx_t id;   /**< The linear id at the level. */
  t8_eclass_t eclass;  /**< The eclass of the tree. */

  /** Order by tree, level and id (the eclass follows from the tree). */
  friend bool
  operator< (const t8_forest_region_element &a, const t8_forest_region_element &b)
  {
    return std::tie (a.gtreeid, a.level, a.id) < std::tie (b.gtreeid, b.level, b.id);
  }
};

/** Insert the descendants of an element of level \a level that touch a boundary of it into the region and the front,
 * recursively over the children touching the boundary.
 * \param [in] scheme   The scheme.
 * \param [in] gtreeid  The global tree of the element.
 * \param [in] eclass   The eclass of the tree.
 * \param [in] element  An element not finer than \a level.
 * \param [in] bdy_dim  The dimension of a boundary of \a element.
 * \param [in] bdy_id   The id of the boundary.
 * \param [in] level    The level of the region elements.
 * \param [in,out] region  The region.
 * \param [in,out] front   The elements inserted into the region are appended.
 */
static void
t8_forest_overlap_boundary_descendants (const t8_scheme *scheme, const t8_gloidx_t gtreeid, const t8_eclass_t eclass,
                                        const t8_element_t *element, const int bdy_dim, const int bdy_id,
                                        const int level, std::set<t8_forest_region_element> &region,
                                        std::vector<t8_forest_region_element> &front)
{
  if (scheme->element_get_level (eclass, element) == level) {
    const t8_forest_region_element region_element
      = { gtreeid, level, scheme->element_get_linear_id (eclass, element, level), eclass };
    if (region.insert (region_element).second) {
      front.push_back (region_element);
    }
    return;
  }
  const int num_children = scheme->element_get_num_boundary_touching_children (eclass, element, bdy_dim, bdy_id);
  std::vector<t8_element_t *> children (num_children);
  std::vector<int> child_bdy_dims (num_children), child_bdy_ids (num_children);
  scheme->element_new (eclass, num_children, children.data ());
  scheme->element_get_children_touching_boundary (eclass, element, bdy_dim, bdy_id, children.data (), num_children,
                                                  child_bdy_dims.data (), child_bdy_ids.data ());
  for (int ichild = 0; ichild < num_children; ichild++) {
    t8_forest_overlap_boundary_descendants (scheme, gtreeid, eclass, children[ichild], child_bdy_dims[ichild],
                                            child_bdy_ids[ichild], level, region, front);
  }
  scheme->element_destroy (eclass, num_children, children.data ());
}

/** The region of this rank: the elements of level \a level that intersect the part of this rank and touch its
 * boundary, grown by \a layers layers of vertex neighbors of the same level. A cover element finer than the level is
 * replaced by its ancestor, a coarser one by its descendants of the level that touch its boundary (the others lie in
 * the part, as do their neighbors).
 * \param [in] forest  A committed forest with splitters.
 * \param [in] level   The level of the region elements.
 * \param [in] layers  The number of layers.
 * \return The elements of the region, sorted by tree, level and id.
 */
static std::set<t8_forest_region_element>
t8_forest_overlap_region (t8_forest_t forest, const int level, const int layers)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  std::set<t8_forest_region_element> region;
  std::vector<t8_forest_region_element> front;

  /* The cover elements at the level: the ancestors of finer ones, the boundary descendants of coarser ones */
  std::vector<std::pair<t8_locidx_t, t8_element_t *>> cover;
  t8_forest_partition_cover (forest, cover);
  for (auto &[ltreeid, element] : cover) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
    const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);
    if (scheme->element_get_level (eclass, element) >= level) {
      t8_element_t *ancestor;
      scheme->element_new (eclass, 1, &ancestor);
      scheme->element_copy (eclass, element, ancestor);
      while (scheme->element_get_level (eclass, ancestor) > level) {
        scheme->element_get_parent (eclass, ancestor, ancestor);
      }
      t8_forest_overlap_boundary_descendants (scheme, gtreeid, eclass, ancestor, 0, 0, level, region, front);
      scheme->element_destroy (eclass, 1, &ancestor);
    }
    else {
      /* Every descendant touching the boundary of the element touches one of its faces */
      const int dim = t8_eclass_to_dimension[eclass];
      const int num_faces = scheme->element_get_num_boundaries (eclass, element, dim - 1);
      for (int face = 0; face < num_faces; face++) {
        t8_forest_overlap_boundary_descendants (scheme, gtreeid, eclass, element, dim - 1, face, level, region,
                                                front);
      }
    }
    scheme->element_destroy (eclass, 1, &element);
  }

  /* The layers: the elements of the same level sharing a vertex (and so also an edge or a face) */
  sc_array_t neighbors;
  sc_array_init (&neighbors, sizeof (t8_forest_boundary_neighbor_t));
  for (int layer = 0; layer < layers; layer++) {
    std::vector<t8_forest_region_element> next;
    for (const t8_forest_region_element &region_element : front) {
      const t8_locidx_t cmesh_ltreeid = t8_cmesh_get_local_id (cmesh, region_element.gtreeid);
      SC_CHECK_ABORT (cmesh_ltreeid >= 0, "The overlap region reaches a tree that the coarse mesh of this rank does "
                                          "not store. Use a replicated coarse mesh or fewer layers.");
      const t8_eclass_t eclass = region_element.eclass;
      t8_element_t *element;
      scheme->element_new (eclass, 1, &element);
      scheme->element_set_linear_id (eclass, element, region_element.level, region_element.id);
      const int num_vertices = scheme->element_get_num_boundaries (eclass, element, 0);
      for (int vertex = 0; vertex < num_vertices; vertex++) {
        t8_forest_element_boundary_neighbors_cmesh (forest, cmesh_ltreeid, region_element.gtreeid, eclass, element, 0,
                                                    vertex, &neighbors);
        for (size_t ineigh = 0; ineigh < neighbors.elem_count; ineigh++) {
          const t8_forest_boundary_neighbor_t *neighbor
            = (const t8_forest_boundary_neighbor_t *) sc_array_index (&neighbors, ineigh);
          const int neigh_level = scheme->element_get_level (neighbor->eclass, neighbor->neigh);
          const t8_forest_region_element neigh_element
            = { neighbor->gtreeid, neigh_level,
                scheme->element_get_linear_id (neighbor->eclass, neighbor->neigh, neigh_level), neighbor->eclass };
          if (region.insert (neigh_element).second) {
            next.push_back (neigh_element);
          }
        }
        t8_forest_boundary_neighbors_clear (forest, &neighbors);
      }
      scheme->element_destroy (eclass, 1, &element);
    }
    front = std::move (next);
  }
  sc_array_reset (&neighbors);
  return region;
}

/** Append the local leaves of a tree that intersect an element: the leaves inside it and the leaf containing it.
 * \param [in] forest     A committed forest.
 * \param [in] ltreeid    A local tree.
 * \param [in] eclass     Its eclass.
 * \param [in] first      The position of the first descendant of the element, in the tree.
 * \param [in] last       The position of the last descendant of the element.
 * \param [in] level      The level of the element.
 * \param [in] id         The linear id of the element at its level.
 * \param [in,out] leaf_ids  For each leaf of the tree, the linear id of its first descendant at the maxlevel, or
 *                           the maximal linear id until it is needed.
 * \param [in,out] leaves    The indices of the leaves in the tree are appended.
 */
static void
t8_forest_leaves_intersecting (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_eclass_t eclass,
                               const t8_forest_position &first, const t8_forest_position &last, const int level,
                               const t8_linearidx_t id, std::vector<t8_linearidx_t> &leaf_ids,
                               std::vector<t8_locidx_t> &leaves)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_element_array_t *tree_leaves = t8_forest_tree_get_leaf_elements (forest, ltreeid);
  const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest, ltreeid);
  const auto leaf_id = [&] (const t8_locidx_t ileaf) {
    if (leaf_ids[ileaf] == std::numeric_limits<t8_linearidx_t>::max ()) {
      leaf_ids[ileaf]
        = scheme->element_get_linear_id (eclass, t8_element_array_index_locidx (tree_leaves, ileaf), forest->maxlevel);
    }
    return leaf_ids[ileaf];
  };
  /* The first leaf whose first descendant is not before id (or after it, if inclusive) */
  const auto search = [&] (const t8_linearidx_t search_id, const bool inclusive) {
    t8_locidx_t lo = 0, hi = num_leaves;
    while (lo < hi) {
      const t8_locidx_t mid = lo + (hi - lo) / 2;
      const t8_linearidx_t mid_id = leaf_id (mid);
      if (mid_id < search_id || (inclusive && mid_id == search_id)) {
        lo = mid + 1;
      }
      else {
        hi = mid;
      }
    }
    return lo;
  };
  const t8_locidx_t begin = search (first.id, false);
  const t8_locidx_t end = search (last.id, true);
  for (t8_locidx_t ileaf = begin; ileaf < end; ileaf++) {
    leaves.push_back (ileaf);
  }
  if (begin == end && begin > 0) {
    /* No leaf starts in the element, the one before may contain it */
    const t8_element_t *leaf = t8_element_array_index_locidx (tree_leaves, begin - 1);
    const int leaf_level = scheme->element_get_level (eclass, leaf);
    if (leaf_level < level) {
      t8_element_t *element;
      scheme->element_new (eclass, 1, &element);
      scheme->element_set_linear_id (eclass, element, level, id);
      if (scheme->element_get_linear_id (eclass, leaf, leaf_level)
          == scheme->element_get_linear_id (eclass, element, leaf_level)) {
        leaves.push_back (begin - 1);
      }
      scheme->element_destroy (eclass, 1, &element);
    }
  }
}

void
t8_forest_ghost_overlap_mirrors (t8_forest_t forest, const int level, const int layers,
                                 std::vector<t8_forest_mirror> &mirrors, std::vector<int> &senders)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  SC_CHECK_ABORT (level >= 0 && layers >= 0, "The level and the number of layers of the overlap must not be negative.");
  mirrors.clear ();
  senders.clear ();
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  if (layers > 0) {
    /* The neighbors across tree boundaries need the vertex and edge connectivity of the trees */
    t8_cmesh_build_boundary_connectivity (t8_forest_get_cmesh (forest), forest->mpicomm);
  }

  /* The elements of the region that do not lie in the own part */
  std::vector<t8_forest_overlap_query> queries;
  if (splitters.has_elements) {
    for (const t8_forest_region_element &region_element : t8_forest_overlap_region (forest, level, layers)) {
      t8_element_t *element;
      scheme->element_new (region_element.eclass, 1, &element);
      scheme->element_set_linear_id (region_element.eclass, element, region_element.level, region_element.id);
      t8_forest_overlap_query query;
      t8_forest_element_positions (forest, region_element.gtreeid, region_element.eclass, element, &query.first,
                                   &query.last);
      scheme->element_destroy (region_element.eclass, 1, &element);
      if (splitters.first_position <= query.first && query.last < splitters.next_position) {
        continue;
      }
      query.id = region_element.id;
      query.level = region_element.level;
      query.source = forest->mpirank;
      query.source_first = splitters.first_position;
      queries.push_back (query);
    }
  }

  /* Route them to the ranks that own a part of them */
  t8_forest_route (forest, queries, [] (const t8_forest_overlap_query &query, const t8_forest_bisection &bisection) {
    return t8_forest_route_halves_positions (bisection, query.first, query.last, false);
  });

  /* The local leaves intersecting the received elements */
  std::vector<std::vector<t8_linearidx_t>> leaf_ids (t8_forest_get_num_local_trees (forest));
  std::vector<t8_locidx_t> leaves;
  std::vector<t8_forest_overlap_notice> notices;
  for (const t8_forest_overlap_query &query : queries) {
    if (query.source == forest->mpirank) {
      continue;
    }
    /* The element lies in one tree; the routing delivers it to the ranks owning a part of it */
    const t8_locidx_t ltreeid = query.first.gtreeid - forest->first_local_tree;
    if (ltreeid < 0 || ltreeid >= t8_forest_get_num_local_trees (forest)) {
      continue;
    }
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltreeid);
    if (leaf_ids[ltreeid].empty ()) {
      leaf_ids[ltreeid].assign (t8_forest_get_tree_num_leaf_elements (forest, ltreeid),
                                std::numeric_limits<t8_linearidx_t>::max ());
    }
    leaves.clear ();
    t8_forest_leaves_intersecting (forest, ltreeid, eclass, query.first, query.last, query.level, query.id,
                                   leaf_ids[ltreeid], leaves);
    for (const t8_locidx_t ileaf : leaves) {
      mirrors.push_back ({ query.source, ltreeid, ileaf });
    }
    if (!leaves.empty ()) {
      notices.push_back ({ query.source_first, forest->mpirank });
    }
  }
  std::sort (mirrors.begin (), mirrors.end ());
  mirrors.erase (std::unique (mirrors.begin (), mirrors.end ()), mirrors.end ());

  /* Tell the receivers of overlap ghosts who sends them, at the owners of their first positions */
  std::sort (notices.begin (), notices.end ());
  notices.erase (std::unique (notices.begin (), notices.end ()), notices.end ());
  t8_forest_route (forest, notices, [] (const t8_forest_overlap_notice &notice, const t8_forest_bisection &bisection) {
    return t8_forest_route_halves_positions (bisection, notice.target, notice.target, true);
  });
  for (const t8_forest_overlap_notice &notice : notices) {
    if (notice.sender != forest->mpirank) {
      senders.push_back (notice.sender);
    }
  }
  std::sort (senders.begin (), senders.end ());
  senders.erase (std::unique (senders.begin (), senders.end ()), senders.end ());
}

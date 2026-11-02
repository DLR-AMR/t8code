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

/** \file t8_forest_balance.cxx
 * Implements functions declared in \ref t8_forest_balance.h.
 */

#include <t8.h>
#include <sc_statistics.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_splitters.hxx>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_cmesh/t8_cmesh.hxx>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/** An element of any tree at a level known from the context: its global tree and its linear id at that level. */
struct t8_balance_key
{
  t8_gloidx_t gtreeid; /**< The global id of the tree. */
  t8_linearidx_t id;   /**< The linear id of the element at its level. */

  /** Equal tree and linear id. */
  bool
  operator== (const t8_balance_key &other) const
  {
    return gtreeid == other.gtreeid && id == other.id;
  }
};

/** A hash of \ref t8_balance_key */
struct t8_balance_key_hash
{
  /** Combine the hashes of the tree and the linear id. */
  size_t
  operator() (const t8_balance_key &key) const
  {
    return std::hash<t8_linearidx_t> () (key.id) ^ (std::hash<t8_gloidx_t> () (key.gtreeid) * 0x9e3779b97f4a7c15ULL);
  }
};

/** For each level, the elements that balance refines, each with whether one of its children is refined as well. */
using t8_balance_refined_t = std::vector<std::unordered_map<t8_balance_key, bool, t8_balance_key_hash>>;

/** For each level, the elements that balance refines. */
using t8_balance_refine_t = std::vector<std::unordered_set<t8_balance_key, t8_balance_key_hash>>;

/** A refined element in the part of another rank, routed to its owner. */
struct t8_balance_item
{
  t8_forest_position first; /**< The position of the first descendant of the element. */
  t8_forest_position last;  /**< The position of the last descendant of the element. */
  t8_linearidx_t id;        /**< The linear id of the element at its level. */
  int level;                /**< The level of the element. */
};

/** The elements that the leaves of this rank force balance to refine, in all trees: the strict ancestors of the
 * leaves, the parents of refined elements, and the elements touching a refined element of the next finer level with
 * a boundary of dimension at least \a min_dim. The levels are processed from the finest to the coarsest.
 *
 * Performance (benchmarks/t8_time_ghost_balance, spherical shell, hexahedra): balance is faster than p4est on 4 and 8
 * ranks for faces, edges and vertices, and 2.5 to 3 times slower on one rank. On one rank, a quarter of the time goes
 * to linear ids of the standalone scheme, which loop over the levels, and a sixth each to the neighbors and to the
 * recursive adapt. Ideas to speed this up:
 * - Inside the own range, the internal nodes of the input forest are refined anyway, so new elements only arise
 *   next to input leaves at least two levels coarser (the seeds of p4est's balance). Find them by locating the leaf
 *   that contains the first descendant of a neighbor (binary search in the leaf array), instead of inserting the
 *   parents of all neighbors of all internal nodes into the hash sets.
 * - Reuse element buffers and a boundary scratch instead of element_new/destroy per element and neighbor, and use
 *   element_get_boundary_neighbors_inside directly for boundaries that do not touch the root boundary.
 * \param [in] forest   A committed forest with a replicated coarse mesh.
 * \param [in] min_dim  The minimal dimension of the boundaries of balance, only faces if it is the dimension minus one.
 * \return For each level, the refined elements, each with whether one of its children is refined.
 */
static t8_balance_refined_t
t8_forest_balance_closure (t8_forest_t forest, const int min_dim)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const int dim = forest->dimension;
  const bool faces_only = min_dim == dim - 1;

  /* The parents of the leaves, their ancestors follow as parents of refined elements */
  t8_balance_refined_t refined (1);
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, itree);
    const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, itree);
    t8_element_t *parent;
    scheme->element_new (eclass, 1, &parent);
    const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest, itree);
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ileaf++) {
      const t8_element_t *leaf = t8_forest_get_leaf_element_in_tree (forest, itree, ileaf);
      const int level = scheme->element_get_level (eclass, leaf);
      if (level == 0) {
        continue;
      }
      if (level > (int) refined.size ()) {
        refined.resize (level);
      }
      scheme->element_get_parent (eclass, leaf, parent);
      refined[level - 1].emplace (t8_balance_key { gtreeid, scheme->element_get_linear_id (eclass, parent, level - 1) },
                                  false);
    }
    scheme->element_destroy (eclass, 1, &parent);
  }

  /* From the finest level to the coarsest, add the parents of the refined elements and the elements touching them.
   * An element of the coarser level touches a refined child of p with dimension at least min_dim if and only if it is
   * p or a neighbor of p at a boundary of p that the child touches with this dimension. So the neighbors are computed
   * once per parent of refined elements, and only at its boundaries touched by a refined child. */
  sc_array_t neighbors;
  sc_array_init (&neighbors, sizeof (t8_forest_boundary_neighbor_t));
  std::vector<t8_balance_key> touching;
  std::vector<t8_element_t *> children;
  std::vector<int> child_bdy_dims;
  for (int level = (int) refined.size () - 1; level >= 1; level--) {
    /* The parents of the refined elements, which have a refined child */
    for (const auto &[key, has_refined_child] : refined[level]) {
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, t8_cmesh_get_local_id (cmesh, key.gtreeid));
      t8_element_t *element, *parent;
      scheme->element_new (eclass, 1, &element);
      scheme->element_new (eclass, 1, &parent);
      scheme->element_set_linear_id (eclass, element, level, key.id);
      scheme->element_get_parent (eclass, element, parent);
      refined[level - 1][{ key.gtreeid, scheme->element_get_linear_id (eclass, parent, level - 1) }] = true;
      scheme->element_destroy (eclass, 1, &element);
      scheme->element_destroy (eclass, 1, &parent);
    }
    /* The neighbors of these parents at the boundaries touched by a refined child */
    touching.clear ();
    for (const auto &[key, has_refined_child] : refined[level - 1]) {
      if (!has_refined_child) {
        continue;
      }
      const t8_locidx_t lctreeid = t8_cmesh_get_local_id (cmesh, key.gtreeid);
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, lctreeid);
      t8_element_t *parent;
      scheme->element_new (eclass, 1, &parent);
      scheme->element_set_linear_id (eclass, parent, level - 1, key.id);
      /* Whether one of the children is refined */
      const auto any_refined = [&] (t8_element_t **children, const int num_children, const int *child_bdy_dims) {
        for (int ichild = 0; ichild < num_children; ichild++) {
          if ((child_bdy_dims == NULL || child_bdy_dims[ichild] >= min_dim)
              && refined[level].count (
                { key.gtreeid, scheme->element_get_linear_id (eclass, children[ichild], level) })) {
            return true;
          }
        }
        return false;
      };
      if (faces_only) {
        const int num_faces = scheme->element_get_num_faces (eclass, parent);
        for (int face = 0; face < num_faces; face++) {
          const int num_children = scheme->element_get_num_face_children (eclass, parent, face);
          children.resize (num_children);
          scheme->element_new (eclass, num_children, children.data ());
          scheme->element_get_children_at_face (eclass, parent, face, children.data (), num_children, NULL);
          const bool touched = any_refined (children.data (), num_children, NULL);
          scheme->element_destroy (eclass, num_children, children.data ());
          if (!touched) {
            continue;
          }
          const t8_eclass_t neigh_eclass
            = t8_forest_element_neighbor_eclass_cmesh (forest, lctreeid, eclass, parent, face);
          t8_element_t *neigh;
          scheme->element_new (neigh_eclass, 1, &neigh);
          int dual_face;
          const t8_gloidx_t neigh_gtreeid = t8_forest_element_face_neighbor_cmesh (
            forest, lctreeid, key.gtreeid, eclass, parent, neigh, neigh_eclass, face, &dual_face);
          if (neigh_gtreeid >= 0) {
            touching.push_back ({ neigh_gtreeid, scheme->element_get_linear_id (neigh_eclass, neigh, level - 1) });
          }
          scheme->element_destroy (neigh_eclass, 1, &neigh);
        }
      }
      else {
        for (int bdy_dim = min_dim; bdy_dim < dim; bdy_dim++) {
          const int num_boundaries = scheme->element_get_num_boundaries (eclass, parent, bdy_dim);
          for (int bdy_id = 0; bdy_id < num_boundaries; bdy_id++) {
            const int num_children
              = scheme->element_get_num_boundary_touching_children (eclass, parent, bdy_dim, bdy_id);
            children.resize (num_children);
            child_bdy_dims.resize (num_children);
            scheme->element_new (eclass, num_children, children.data ());
            scheme->element_get_children_touching_boundary (eclass, parent, bdy_dim, bdy_id, children.data (),
                                                            num_children, child_bdy_dims.data (), NULL);
            const bool touched = any_refined (children.data (), num_children, child_bdy_dims.data ());
            scheme->element_destroy (eclass, num_children, children.data ());
            if (!touched) {
              continue;
            }
            t8_forest_element_boundary_neighbors_cmesh (forest, lctreeid, key.gtreeid, eclass, parent, bdy_dim, bdy_id,
                                                        &neighbors);
            for (size_t ineigh = 0; ineigh < neighbors.elem_count; ineigh++) {
              const t8_forest_boundary_neighbor_t *neighbor
                = (const t8_forest_boundary_neighbor_t *) sc_array_index (&neighbors, ineigh);
              touching.push_back (
                { neighbor->gtreeid, scheme->element_get_linear_id (neighbor->eclass, neighbor->neigh, level - 1) });
            }
            t8_forest_boundary_neighbors_clear (forest, &neighbors);
          }
        }
      }
      scheme->element_destroy (eclass, 1, &parent);
    }
    for (const t8_balance_key &key : touching) {
      refined[level - 1].emplace (key, false);
    }
  }
  sc_array_reset (&neighbors);
  return refined;
}

/** Keep only the routed elements without a descendant among them, sorted by position.
 * \param [in,out] items  The routed elements.
 */
static void
t8_forest_balance_merge_items (std::vector<t8_balance_item> &items)
{
  /* Sort by position, the coarser element first at equal positions, and remove duplicates */
  std::sort (items.begin (), items.end (), [] (const t8_balance_item &a, const t8_balance_item &b) {
    return a.first < b.first || (a.first == b.first && a.level < b.level);
  });
  items.erase (std::unique (items.begin (), items.end (),
                            [] (const t8_balance_item &a, const t8_balance_item &b) {
                              return a.first == b.first && a.level == b.level;
                            }),
               items.end ());
  /* The descendants of an element follow it directly, so it has one if the next element starts inside it */
  std::vector<bool> keep (items.size ());
  for (size_t i = 0; i < items.size (); i++) {
    keep[i] = !(i + 1 < items.size () && items[i + 1].first <= items[i].last);
  }
  size_t num_kept = 0;
  for (size_t i = 0; i < items.size (); i++) {
    if (keep[i]) {
      items[num_kept++] = items[i];
    }
  }
  items.resize (num_kept);
}

/** Adapt callback of the one pass balance: refine the elements in the set of the forest's t8code_data. */
static int
t8_forest_balance_one_pass_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id,
                                  const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id,
                                  const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                                  [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const t8_balance_refine_t *refine = (const t8_balance_refine_t *) forest->t8code_data;
  const int level = scheme->element_get_level (tree_class, elements[0]);
  if (level >= (int) refine->size ()) {
    return 0;
  }
  const t8_balance_key key = { t8_forest_global_tree_id (forest_from, ltree_id),
                               scheme->element_get_linear_id (tree_class, elements[0], level) };
  return (*refine)[level].count (key) > 0;
}

#if T8_ENABLE_DEBUG
/** Whether two forests with the same partition have the same leaves. */
static bool
t8_forest_balance_same_leaves (t8_forest_t forest_a, t8_forest_t forest_b)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest_a);
  if (t8_forest_get_num_local_trees (forest_a) != t8_forest_get_num_local_trees (forest_b)) {
    t8_errorf ("Balance differs in the number of trees: %i, ripple %i\n", t8_forest_get_num_local_trees (forest_a),
               t8_forest_get_num_local_trees (forest_b));
    return false;
  }
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest_a); itree++) {
    const t8_locidx_t num_leaves = t8_forest_get_tree_num_leaf_elements (forest_a, itree);
    if (num_leaves != t8_forest_get_tree_num_leaf_elements (forest_b, itree)) {
      t8_errorf ("Balance differs in the number of leaves of tree %lli: %i, ripple %i\n",
                 (long long) t8_forest_global_tree_id (forest_a, itree), num_leaves,
                 t8_forest_get_tree_num_leaf_elements (forest_b, itree));
      const t8_eclass_t eclass = t8_forest_get_tree_class (forest_a, itree);
      const t8_locidx_t num_common = SC_MIN (num_leaves, t8_forest_get_tree_num_leaf_elements (forest_b, itree));
      for (t8_locidx_t ileaf = 0; ileaf < num_common; ileaf++) {
        const t8_element_t *leaf_a = t8_forest_get_leaf_element_in_tree (forest_a, itree, ileaf);
        const t8_element_t *leaf_b = t8_forest_get_leaf_element_in_tree (forest_b, itree, ileaf);
        if (!scheme->element_is_equal (eclass, leaf_a, leaf_b)) {
          const int level_a = scheme->element_get_level (eclass, leaf_a);
          const int level_b = scheme->element_get_level (eclass, leaf_b);
          t8_errorf ("First difference at leaf %i: level %i id %llu, ripple level %i id %llu\n", ileaf, level_a,
                     (unsigned long long) scheme->element_get_linear_id (eclass, leaf_a, level_a), level_b,
                     (unsigned long long) scheme->element_get_linear_id (eclass, leaf_b, level_b));
          break;
        }
      }
      return false;
    }
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest_a, itree);
    for (t8_locidx_t ileaf = 0; ileaf < num_leaves; ileaf++) {
      const t8_element_t *leaf_a = t8_forest_get_leaf_element_in_tree (forest_a, itree, ileaf);
      const t8_element_t *leaf_b = t8_forest_get_leaf_element_in_tree (forest_b, itree, ileaf);
      if (!scheme->element_is_equal (eclass, leaf_a, leaf_b)) {
        const int level_a = scheme->element_get_level (eclass, leaf_a);
        const int level_b = scheme->element_get_level (eclass, leaf_b);
        t8_errorf ("Balance differs in tree %lli at leaf %i: level %i id %llu, ripple level %i id %llu\n",
                   (long long) t8_forest_global_tree_id (forest_a, itree), ileaf, level_a,
                   (unsigned long long) scheme->element_get_linear_id (eclass, leaf_a, level_a), level_b,
                   (unsigned long long) scheme->element_get_linear_id (eclass, leaf_b, level_b));
        return false;
      }
    }
  }
  return true;
}
#endif

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/**
 * Check whether the leaves at a face of an element are more than one level finer than the element, i.e. whether
 * one of its half face neighbors has local or ghost leaf descendants. Works across tree boundaries.
 *
 * \param[in] forest_from  The forest to check.
 * \param[in] ltree_id     The local id of the tree of \a element.
 * \param[in] tree_class   The class of the tree.
 * \param[in] scheme       The scheme of the forest.
 * \param[in] element      The element.
 * \param[in] iface        A face of \a element.
 * \return 1 if \a element has to be refined because of the face, 0 otherwise.
 */
static int
t8_forest_balance_face_is_unbalanced (t8_forest_t forest_from, const t8_locidx_t ltree_id, const t8_eclass_t tree_class,
                                      const t8_scheme *scheme, const t8_element_t *element, const int iface)
{
  int unbalanced = 0;
  /* Get the element class and scheme of the face neighbor */
  const t8_eclass_t neigh_class = t8_forest_element_neighbor_eclass (forest_from, ltree_id, element, iface);
  /* Allocate memory for the number of half face neighbors */
  const int num_half_neighbors = scheme->element_get_num_face_children (tree_class, element, iface);
  t8_element_t **half_neighbors = T8_ALLOC (t8_element_t *, num_half_neighbors);
  scheme->element_new (neigh_class, num_half_neighbors, half_neighbors);
  /* Compute the half face neighbors of element at this face */
  const t8_gloidx_t neighbor_tree = t8_forest_element_half_face_neighbors (
    forest_from, ltree_id, element, half_neighbors, neigh_class, iface, num_half_neighbors, NULL);
  if (neighbor_tree >= 0) {
    /* The face neighbors do exist, check for each one, whether it has
     * local or ghost leaf descendants in the forest.
     * If so, the element will be refined. */
    for (int ineigh = 0; ineigh < num_half_neighbors && !unbalanced; ineigh++) {
      unbalanced = t8_forest_element_has_leaf_desc (forest_from, neighbor_tree, half_neighbors[ineigh], neigh_class);
    }
  }
  /* clean-up */
  scheme->element_destroy (neigh_class, num_half_neighbors, half_neighbors);
  T8_FREE (half_neighbors);
  return unbalanced;
}

/** Scratch arrays for t8_forest_balance_adapt, set during a balance round and in t8_forest_is_balanced, so that
 * the neighbor arrays are not allocated for every boundary of every element. NULL otherwise. */
static t8_forest_boundary_scratch_t *balance_scratch = NULL;

/**
 * Check whether leaves that touch a boundary of an element with a boundary of dimension at least min_dim are more
 * than one level finer than the element, i.e. whether a child of a neighbor (of the same level) at this boundary,
 * which touches it with dimension at least min_dim, has local or ghost leaf descendants. Then all leaves of this
 * child, including the ones touching the boundary, are at least two levels finer. The child may touch the boundary
 * in a part that is not a boundary of the element, e.g. at the midpoint of an edge of a tetrahedron.
 * The neighbors are those of \ref t8_forest_element_boundary_neighbors, inside the tree and in the trees at the
 * root boundary that contains the boundary.
 *
 * \param[in] forest_from  The forest to check.
 * \param[in] ltree_id     The local id of the tree of \a element.
 * \param[in] tree_class   The class of the tree.
 * \param[in] scheme       The scheme of the forest.
 * \param[in] element      The element.
 * \param[in] bdy_dim      The dimension of the boundary.
 * \param[in] bdy_id       The boundary of \a element of dimension \a bdy_dim.
 * \param[in] min_dim      The minimum dimension of the intersection of a leaf with the boundary.
 * \param[in,out] neighbors An initialized array of \ref t8_forest_boundary_neighbor_t for the neighbors (empty on
 *                          input and output).
 * \return 1 if \a element has to be refined because of the boundary, 0 otherwise.
 */
static int
t8_forest_balance_boundary_is_unbalanced (t8_forest_t forest_from, const t8_locidx_t ltree_id,
                                          [[maybe_unused]] const t8_eclass_t tree_class, const t8_scheme *scheme,
                                          const t8_element_t *element, const int bdy_dim, const int bdy_id,
                                          const int min_dim, sc_array_t *neighbors)
{
  int unbalanced = 0;
  t8_forest_element_boundary_neighbors (forest_from, ltree_id, element, bdy_dim, bdy_id, neighbors);
  for (size_t ineigh = 0; ineigh < neighbors->elem_count && !unbalanced; ineigh++) {
    const t8_forest_boundary_neighbor_t *neighbor
      = (const t8_forest_boundary_neighbor_t *) sc_array_index (neighbors, ineigh);
    const t8_eclass_t neigh_class = neighbor->eclass;
    /* If the neighbor has no finer leaves, none of its children has, so we skip them.
     * This is the common case and much cheaper than checking each child. */
    if (!t8_forest_element_has_leaf_desc (forest_from, neighbor->gtreeid, neighbor->neigh, neigh_class)) {
      continue;
    }
    const int max_children = scheme->get_max_num_boundary_touching_children (neigh_class);
    std::vector<t8_element_t *> children (max_children);
    std::vector<int> child_bdy_dims (max_children);
    scheme->element_new (neigh_class, max_children, children.data ());
    const int num_children
      = scheme->element_get_num_boundary_touching_children (neigh_class, neighbor->neigh, bdy_dim, neighbor->bdy_id);
    scheme->element_get_children_touching_boundary (neigh_class, neighbor->neigh, bdy_dim, neighbor->bdy_id,
                                                    children.data (), num_children, child_bdy_dims.data (), NULL);
    for (int ichild = 0; ichild < num_children && !unbalanced; ichild++) {
      if (child_bdy_dims[ichild] >= min_dim) {
        unbalanced = t8_forest_element_has_leaf_desc (forest_from, neighbor->gtreeid, children[ichild], neigh_class);
      }
    }
    scheme->element_destroy (neigh_class, max_children, children.data ());
  }
  t8_forest_boundary_neighbors_clear (forest_from, neighbors);
  return unbalanced;
}

/** This is the adapt function called during one round of balance.
 * We refine an element if it has any face neighbor with a level larger
 * than the element's level + 1.
 *
 * TODO: We currently do not adapt recursively since some functions such
 * as half neighbor computation require the forest to be committed. Thus,
 * we pass forest_from as a parameter. But doing so is not valid anymore
 * if we refine recursively.
 *
 * \param[in, out] forest       The forest to be adapted / balanced.
 * \param[in]      forest_from  The forest from which the current one is derived.
 * \param[in]      ltree_id     The local id of the tree the element is in.
 * \param[in]      tree_class   The element class of the tree the element is in.
 * \param[in]      lelement_id  The local id of the element within the tree.
 * \param[in]      scheme       The scheme class.
 * \param[in]      is_family    A switch indicating whether the passed elements form a family.
 * \param[in]      num_elements The number of elements passed as input.
 * \param[in]      elements     The elements array.
 *
 * \return 1 if the element(s) has/have to be refined, 0 otherwise.
 */
static int
t8_forest_balance_adapt (t8_forest_t forest, t8_forest_t forest_from, const t8_locidx_t ltree_id,
                         const t8_eclass_t tree_class, [[maybe_unused]] const t8_locidx_t lelement_id,
                         const t8_scheme *scheme, [[maybe_unused]] const int is_family,
                         [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
  const t8_element_t *element = elements[0];

  /* We only need to check an element, if its level is smaller then the maximum
   * level in the forest minus 2.
   * Otherwise there cannot exist neighbors of level greater than the element's level plus one.
   * The variable maxlevel_existing is only set if we enter this function from t8_forest_balance.
   * If we enter from the check function is_balanced, then it may not be set.
   */

  if (forest_from->maxlevel_existing <= 0
      || scheme->element_get_level (tree_class, element) <= forest_from->maxlevel_existing - 2) {

    int *pdone = (int *) forest->t8code_data;

    /* The smallest dimension of the boundaries considered by the balance type */
    const int dim = t8_eclass_to_dimension[tree_class];
    const int min_bdy_dim = forest->balance_type == T8_GHOST_FACES   ? SC_MAX (dim - 1, 0)
                            : forest->balance_type == T8_GHOST_EDGES ? SC_MAX (dim - 2, 0)
                                                                     : 0;
    if (min_bdy_dim == dim - 1) {
      const int num_faces = scheme->element_get_num_faces (tree_class, element);
      for (int iface = 0; iface < num_faces; iface++) {
        if (t8_forest_balance_face_is_unbalanced (forest_from, ltree_id, tree_class, scheme, element, iface)) {
          /* This element should be refined */
          *pdone = 0;
          return 1;
        }
      }
    }
    else {
      t8_forest_boundary_scratch_t local_scratch;
      if (balance_scratch == NULL) {
        t8_forest_boundary_scratch_init (&local_scratch);
      }
      sc_array_t *neighbors = balance_scratch != NULL ? &balance_scratch->neighbors : &local_scratch.neighbors;
      int unbalanced = 0;
      if (min_bdy_dim < dim - 1) {
        /* Every neighbor (of the same level) at a boundary of the element is also a
         * neighbor at the boundaries of dimension min_bdy_dim contained in this boundary
         * (e.g. a face neighbor at the edges of the face). If no neighbor at a boundary of
         * dimension min_bdy_dim has finer leaves, no boundary is unbalanced, which is the
         * common case and much cheaper than computing the neighbors at all boundaries. */
        bool any_finer = false;
        const int num_min_boundaries = scheme->element_get_num_boundaries (tree_class, element, min_bdy_dim);
        for (int ibdy = 0; ibdy < num_min_boundaries && !any_finer; ibdy++) {
          t8_forest_element_boundary_neighbors (forest_from, ltree_id, element, min_bdy_dim, ibdy, neighbors);
          for (size_t ineigh = 0; ineigh < neighbors->elem_count && !any_finer; ineigh++) {
            const t8_forest_boundary_neighbor_t *neighbor
              = (const t8_forest_boundary_neighbor_t *) sc_array_index (neighbors, ineigh);
            any_finer
              = t8_forest_element_has_leaf_desc (forest_from, neighbor->gtreeid, neighbor->neigh, neighbor->eclass);
          }
          t8_forest_boundary_neighbors_clear (forest_from, neighbors);
        }
        if (!any_finer) {
          if (balance_scratch == NULL) {
            t8_forest_boundary_scratch_reset (&local_scratch);
          }
          return 0;
        }
      }
      /* Loop over the boundaries of all dimensions considered by the balance type */
      for (int bdy_dim = min_bdy_dim; bdy_dim < dim && !unbalanced; bdy_dim++) {
        const int num_boundaries = scheme->element_get_num_boundaries (tree_class, element, bdy_dim);
        for (int bdy_id = 0; bdy_id < num_boundaries && !unbalanced; bdy_id++) {
          unbalanced = t8_forest_balance_boundary_is_unbalanced (forest_from, ltree_id, tree_class, scheme, element,
                                                                 bdy_dim, bdy_id, min_bdy_dim, neighbors);
        }
      }
      if (balance_scratch == NULL) {
        t8_forest_boundary_scratch_reset (&local_scratch);
      }
      if (unbalanced) {
        /* This element should be refined */
        *pdone = 0;
        return 1;
      }
    }
  }

  return 0;
}

/**
 * Collective function to compute the maximum occurring refinement level in a forest
 *
 * \param[in,out] forest  The forest which the maximum refinement level is computed for and stored in.
 */
static void
t8_forest_compute_max_element_level (t8_forest_t forest)
{
  t8_locidx_t ielement, elem_in_tree;
  t8_locidx_t itree, num_trees;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  int local_max_level = 0;

  /* Iterate over all local trees and all local elements and compute the maximum occurring level */
  num_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0; itree < num_trees; itree++) {
    elem_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    for (ielement = 0; ielement < elem_in_tree; ielement++) {
      /* Get the element and compute its level */
      const t8_element_t *elem = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      const int elem_level = scheme->element_get_level (tree_class, elem);
      local_max_level = SC_MAX (local_max_level, elem_level);
    }
  }
  /* Communicate the local maximum levels */
  sc_MPI_Allreduce (&local_max_level, &forest->maxlevel_existing, 1, sc_MPI_INT, sc_MPI_MAX, forest->mpicomm);
}

/** Balance a forest with \ref T8_BALANCE_RIPPLE. */
static void
t8_forest_balance_ripple (t8_forest_t forest, int repartition)
{
  t8_forest_t forest_temp, forest_from, forest_partition;
  int done = 0, done_global = 0;
  int count_rounds = 0;
  /* The following variables are only required if profiling is
   * enabled. */
  int num_stats_allocated, istats;
  const int stat_alloc_chunk_size = 10; /* How many stats we allocate (and add if we need more) */
  int count_adapt_stats = 0, count_ghost_stats = 0;
  int count_partition_stats = 0;
  double ada_time, ghost_time, part_time;
  sc_statinfo_t *adap_stats, *ghost_stats, *partition_stats;

  t8_global_productionf ("Into t8_forest_balance with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_leaf_elements (forest->set_from));
  t8_log_indent_push ();

  /* Set default value to prevent compiler warning */
  adap_stats = ghost_stats = partition_stats = NULL;

  if (forest->profile != NULL) {
    /* Profiling is enable, so we measure the runtime of balance */
    forest->profile->balance_runtime = -sc_MPI_Wtime ();
    /* We store the individual adapt, ghost, and partition runtimes */
    /* We reserve memory for stat_alloc_chunk_size - 1 many balance rounds
     * (the extra entry is required for the total sum).
     * We will grow the statistics arrays dynamically if more rounds are required. */
    T8_ASSERT (stat_alloc_chunk_size > 0);
    num_stats_allocated = stat_alloc_chunk_size;
    adap_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    ghost_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    if (repartition) {
      partition_stats = T8_ALLOC_ZERO (sc_statinfo_t, num_stats_allocated);
    }
  }

  /* Compute the maximum occurring refinement level in the forest */
  t8_forest_compute_max_element_level (forest->set_from);
  t8_global_productionf ("Computed maximum occurring level:\t%i\n", forest->set_from->maxlevel_existing);
  /* Use set_from as the first forest to adapt */
  forest_from = forest->set_from;
  /* This function is reference neutral regarding forest_from */
  t8_forest_ref (forest_from);

  /* The ghost layer must contain all neighbors of the balance type */
  if (forest->set_from->ghosts == NULL || forest->set_from->ghost_type < forest->balance_type) {
    if (forest->set_from->ghosts != NULL) {
      t8_forest_ghost_unref (&forest->set_from->ghosts);
    }
    forest->set_from->ghost_type = forest->balance_type;
    /* The same algorithm as for the ghost layers of the later rounds (the default of t8_forest_set_ghost) */
    t8_forest_ghost_create_touching (forest->set_from);
  }

  while (!done_global) {
    done = 1;

    T8_ASSERT (forest_from->maxlevel_existing >= 0);
    /* Initialize the temp forest to be adapted from forest_from */
    t8_forest_init (&forest_temp);
    /* Update the maximum occurring level */
    forest_temp->maxlevel_existing = forest_from->maxlevel_existing;
    forest_temp->balance_type = forest->balance_type;
    /* Adapt the forest */
    t8_forest_set_adapt (forest_temp, forest_from, t8_forest_balance_adapt, 0);
    if (forest->balance_type != T8_GHOST_FACES) {
      /* The neighbors across tree boundaries need the vertex and edge connectivity of the trees */
      t8_cmesh_build_boundary_connectivity (t8_forest_get_cmesh (forest_from), forest->mpicomm);
    }
    if (!repartition) {
      t8_forest_set_ghost (forest_temp, 1, forest->balance_type);
    }
    forest_temp->t8code_data = &done;
    /* If profiling is enabled, measure ghost/adapt runtimes */
    if (forest->profile != NULL) {
      t8_forest_set_profiling (forest_temp, 1);
    }
    t8_global_productionf ("Profiling: %i\n", forest->profile != NULL);
    /* Adapt the forest */
    /* The adapt function asks whether the same elements (the neighbors of
     * neighboring elements) have leaf descendants many times */
    t8_forest_boundary_scratch_t scratch;
    t8_forest_boundary_scratch_init (&scratch);
    balance_scratch = &scratch;
    t8_forest_leaf_desc_cache_begin (forest_from);
    t8_forest_commit (forest_temp);
    t8_forest_leaf_desc_cache_end ();
    balance_scratch = NULL;
    t8_forest_boundary_scratch_reset (&scratch);
    /* Store the runtimes of adapt and ghost */
    if (forest->profile != NULL) {
      if (count_rounds > num_stats_allocated - 2) {
        T8_ASSERT (count_adapt_stats <= count_rounds);
        T8_ASSERT (count_ghost_stats <= count_rounds);
        T8_ASSERT (count_partition_stats <= count_rounds);
        /* re-allocate memory for stats */
        num_stats_allocated += stat_alloc_chunk_size;
        adap_stats = T8_REALLOC (adap_stats, sc_statinfo_t, num_stats_allocated);
        ghost_stats = T8_REALLOC (ghost_stats, sc_statinfo_t, num_stats_allocated);
        if (repartition) {
          partition_stats = T8_REALLOC (partition_stats, sc_statinfo_t, num_stats_allocated);
        }
      }
      sc_stats_set1 (&adap_stats[count_adapt_stats], forest_temp->profile->adapt_runtime, "forest balance: Adapt time");
      count_adapt_stats++;
      if (!repartition) {
        sc_stats_set1 (&ghost_stats[count_ghost_stats], forest_temp->profile->ghost_runtime,
                       "forest balance: Ghost time");
        count_ghost_stats++;
      }
    }

    /* Compute the logical and of all process local done values, if this results
     * in 1 then all processes are finished */
    sc_MPI_Allreduce (&done, &done_global, 1, sc_MPI_INT, sc_MPI_LAND, forest->mpicomm);

    if (repartition && !done_global) {
      /* If repartitioning is used, we partition the forest */
      t8_forest_init (&forest_partition);
      /* Update the maximum occurring level */
      forest_partition->maxlevel_existing = forest_temp->maxlevel_existing;
      forest_partition->balance_type = forest->balance_type;
      t8_forest_set_partition (forest_partition, forest_temp, 0);
      t8_forest_set_ghost (forest_partition, 1, forest->balance_type);
      /* If profiling is enabled, measure partition runtimes */
      if (forest->profile != NULL) {
        t8_forest_set_profiling (forest_partition, 1);
      }
      t8_forest_commit (forest_partition);

      /* Store the runtimes of partition */
      if (forest->profile != NULL) {
        sc_stats_set1 (&partition_stats[count_partition_stats], forest_partition->profile->partition_runtime,
                       "forest balance: Partition time");
        count_partition_stats++;
        sc_stats_set1 (&ghost_stats[count_ghost_stats], forest_partition->profile->ghost_runtime,
                       "forest balance: Ghost time");
        count_ghost_stats++;
      }

      forest_temp = forest_partition;
      forest_partition = NULL;
    }
    /* Adapt forest_temp in the next round */
    forest_from = forest_temp;
    count_rounds++;
  }

  T8_ASSERT (t8_forest_is_balanced (forest_temp));
  /* Forest_temp is now balanced, we copy its trees and elements to forest */
  t8_forest_copy_trees (forest, forest_temp, 1);
  if (forest->cmesh != forest_temp->cmesh) {
    /* Repartitioning the forest also repartitioned a partitioned cmesh, which has to match the trees */
    t8_cmesh_ref (forest_temp->cmesh);
    t8_cmesh_unref (&forest->cmesh);
    forest->cmesh = forest_temp->cmesh;
  }
  /* TODO: Also copy ghost elements if ghost creation is set */

  t8_log_indent_pop ();
  t8_global_productionf ("Done t8_forest_balance with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_leaf_elements (forest_temp));
  t8_debugf ("t8_forest_balance needed %i rounds.\n", count_rounds);
  /* clean-up */
  t8_forest_unref (&forest_temp);

  if (forest->profile != NULL) {
    /* Profiling is enabled, so we measure the runtime of balance. */
    forest->profile->balance_runtime += sc_MPI_Wtime ();
    forest->profile->balance_rounds = count_rounds;
    /* Print the runtime of adapt/ghost/partition */
    /* Compute the overall runtime and store in last entry */
    ada_time = ghost_time = part_time = 0;
    for (istats = 0; istats < count_adapt_stats; istats++) {
      ada_time += adap_stats[istats].sum_values;
    }
    for (istats = 0; istats < count_ghost_stats; istats++) {
      ghost_time += ghost_stats[istats].sum_values;
    }
    if (repartition) {
      for (istats = 0; istats < count_partition_stats; istats++) {
        part_time += partition_stats[istats].sum_values;
      }
    }
    sc_stats_set1 (&adap_stats[count_adapt_stats], ada_time, "forest balance: Total adapt time");
    sc_stats_set1 (&ghost_stats[count_ghost_stats], ghost_time, "forest balance: Total ghost time");
    if (repartition) {
      sc_stats_set1 (&partition_stats[count_partition_stats], part_time, "forest balance: Total partition time");
    }

    /* Compute and print the intermediate stats */
    T8_ASSERT (count_rounds + 1 <= num_stats_allocated);
    sc_stats_compute (forest->mpicomm, count_adapt_stats + 1, adap_stats);
    sc_stats_compute (forest->mpicomm, count_ghost_stats + 1, ghost_stats);
    if (repartition) {
      sc_stats_compute (forest->mpicomm, count_partition_stats + 1, partition_stats);
    }
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_adapt_stats + 1, adap_stats, 1, 1);
    sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_ghost_stats + 1, ghost_stats, 1, 1);
    if (repartition) {
      sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, count_partition_stats + 1, partition_stats, 1, 1);
    }
    T8_FREE (adap_stats);
    T8_FREE (ghost_stats);
    if (repartition) {
      T8_FREE (partition_stats);
    }
  }
}

/** Balance a forest with \ref T8_BALANCE_ONE_PASS. */
static void
t8_forest_balance_one_pass (t8_forest_t forest, int repartition)
{
  t8_forest_t forest_from = forest->set_from;
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest_from);
  SC_CHECK_ABORT (!t8_cmesh_is_partitioned (cmesh), "Balance in one pass needs a replicated coarse mesh.\n");
  t8_global_productionf ("Into t8_forest_balance (one pass) with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_leaf_elements (forest_from));
  t8_log_indent_push ();
  if (forest->profile != NULL) {
    forest->profile->balance_runtime = -sc_MPI_Wtime ();
  }

  /* The elements that the leaves of this rank force to be refined, in all trees */
  const int dim = forest_from->dimension;
  const int min_dim = forest->balance_type == T8_GHOST_FACES   ? SC_MAX (dim - 1, 0)
                      : forest->balance_type == T8_GHOST_EDGES ? SC_MAX (dim - 2, 0)
                                                               : 0;
  if (min_dim < dim - 1) {
    /* The neighbors across tree boundaries need the vertex and edge connectivity of the trees */
    t8_cmesh_build_boundary_connectivity (cmesh, forest_from->mpicomm);
  }
  const t8_balance_refined_t refined = t8_forest_balance_closure (forest_from, min_dim);

  /* Those in the part of this rank are refined here. Of those in the parts of other ranks, the ones without a refined
   * child are sent to their owners. Elements in several parts contain leaves of several ranks and are refined anyway. */
  const t8_forest_splitters &splitters = t8_forest_get_splitters (forest_from);
  const t8_scheme *scheme = t8_forest_get_scheme (forest_from);
  t8_balance_refine_t refine (refined.size ());
  std::vector<t8_balance_item> items;
  for (int level = 0; level < (int) refined.size (); level++) {
    for (const auto &[key, has_refined_child] : refined[level]) {
      const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, t8_cmesh_get_local_id (cmesh, key.gtreeid));
      t8_element_t *element;
      scheme->element_new (eclass, 1, &element);
      scheme->element_set_linear_id (eclass, element, level, key.id);
      t8_balance_item item;
      t8_forest_element_positions (forest_from, key.gtreeid, eclass, element, &item.first, &item.last);
      scheme->element_destroy (eclass, 1, &element);
      if (splitters.has_elements && splitters.first_position <= item.first && item.last < splitters.next_position) {
        refine[level].insert (key);
      }
      else if (!has_refined_child
               && (!splitters.has_elements || item.last < splitters.first_position
                   || !(item.first < splitters.next_position))) {
        item.id = key.id;
        item.level = level;
        items.push_back (item);
      }
    }
  }
  t8_forest_route (
    forest_from, items,
    [] (const t8_balance_item &item, const t8_forest_bisection &level) {
      return t8_forest_route_halves_positions (level, item.first, item.last, true);
    },
    t8_forest_balance_merge_items);

  /* The received elements and their ancestors are refined as well */
  for (const t8_balance_item &item : items) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, t8_cmesh_get_local_id (cmesh, item.first.gtreeid));
    if (item.level >= (int) refine.size ()) {
      refine.resize (item.level + 1);
    }
    t8_element_t *element, *parent;
    scheme->element_new (eclass, 1, &element);
    scheme->element_new (eclass, 1, &parent);
    scheme->element_set_linear_id (eclass, element, item.level, item.id);
    for (int level = item.level; level >= 0; level--) {
      if (!refine[level]
             .insert ({ item.first.gtreeid, scheme->element_get_linear_id (eclass, element, level) })
             .second) {
        break;
      }
      if (level > 0) {
        scheme->element_get_parent (eclass, element, parent);
        scheme->element_copy (eclass, parent, element);
      }
    }
    scheme->element_destroy (eclass, 1, &element);
    scheme->element_destroy (eclass, 1, &parent);
  }

  /* Refine the leaves recursively */
  t8_forest_t forest_temp;
  t8_forest_init (&forest_temp);
  t8_forest_ref (forest_from);
  t8_forest_set_adapt (forest_temp, forest_from, t8_forest_balance_one_pass_adapt, 1);
  forest_temp->t8code_data = &refine;
  forest_temp->balance_type = forest->balance_type;
  if (forest->profile != NULL) {
    t8_forest_set_profiling (forest_temp, 1);
  }
  t8_forest_commit (forest_temp);
#if T8_ENABLE_DEBUG
  {
    /* The ripple algorithm gives the same leaves */
    t8_forest_t forest_ripple;
    t8_forest_init (&forest_ripple);
    t8_forest_set_balance_type (forest_ripple, forest->balance_type);
    t8_forest_set_balance_algorithm (forest_ripple, T8_BALANCE_RIPPLE);
    t8_forest_ref (forest_from);
    t8_forest_set_balance (forest_ripple, forest_from, 1);
    t8_forest_commit (forest_ripple);
    T8_ASSERT (t8_forest_balance_same_leaves (forest_temp, forest_ripple));
    t8_forest_unref (&forest_ripple);
  }
#endif

  if (repartition) {
    t8_forest_t forest_partition;
    t8_forest_init (&forest_partition);
    t8_forest_set_partition (forest_partition, forest_temp, 0);
    if (forest->profile != NULL) {
      t8_forest_set_profiling (forest_partition, 1);
    }
    t8_forest_commit (forest_partition);
    forest_temp = forest_partition;
  }
  /* Copy the trees and elements to forest */
  t8_forest_copy_trees (forest, forest_temp, 1);
  if (forest->cmesh != forest_temp->cmesh) {
    /* Repartitioning the forest also repartitioned a partitioned cmesh, which has to match the trees */
    t8_cmesh_ref (forest_temp->cmesh);
    t8_cmesh_unref (&forest->cmesh);
    forest->cmesh = forest_temp->cmesh;
  }
  t8_log_indent_pop ();
  t8_global_productionf ("Done t8_forest_balance (one pass) with %lli global elements.\n",
                         (long long) t8_forest_get_global_num_leaf_elements (forest_temp));
  t8_forest_unref (&forest_temp);
  if (forest->profile != NULL) {
    forest->profile->balance_runtime += sc_MPI_Wtime ();
    forest->profile->balance_rounds = 1;
  }
}

void
t8_forest_balance (t8_forest_t forest, int repartition)
{
  t8_forest_balance_algorithm_t algorithm = (t8_forest_balance_algorithm_t) forest->balance_algorithm;
  if (algorithm == T8_BALANCE_DEFAULT) {
    algorithm
      = t8_cmesh_is_partitioned (t8_forest_get_cmesh (forest->set_from)) ? T8_BALANCE_RIPPLE : T8_BALANCE_ONE_PASS;
  }
  if (algorithm == T8_BALANCE_ONE_PASS) {
    t8_forest_balance_one_pass (forest, repartition);
  }
  else {
    t8_forest_balance_ripple (forest, repartition);
  }
}

/* Check whether the local elements of a forest are balanced. */
int
t8_forest_is_balanced (t8_forest_t forest)
{
  t8_forest_t forest_from;
  t8_locidx_t num_trees, num_elements;
  t8_locidx_t itree, ielem;
  void *data_temp;
  int dummy_int;

  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_scheme *scheme = t8_forest_get_scheme (forest);

  /* temporarily save forest_from */
  forest_from = forest->set_from;

  forest->set_from = forest;

  /* temporarily save forest t8code_data */
  data_temp = forest->t8code_data;
  forest->t8code_data = &dummy_int;

  /* The scratch arrays for t8_forest_balance_adapt (unless we are inside a balance round) */
  t8_forest_boundary_scratch_t scratch;
  const bool own_scratch = balance_scratch == NULL;
  if (own_scratch) {
    t8_forest_boundary_scratch_init (&scratch);
    balance_scratch = &scratch;
  }

  num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over all trees */
  for (itree = 0; itree < num_trees; itree++) {
    num_elements = t8_forest_get_tree_num_leaf_elements (forest, itree);
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Iterate over all elements of this tree */
    for (ielem = 0; ielem < num_elements; ielem++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielem);
      /* Test if this element would need to be refined in the balance step.
       * If so, the forest is not balanced locally. */
      if (t8_forest_balance_adapt (forest, forest, itree, tree_class, ielem, scheme, 0, 1,
                                   (t8_element_t **) (&element))) {
        forest->set_from = forest_from;
        forest->t8code_data = data_temp;
        if (own_scratch) {
          balance_scratch = NULL;
          t8_forest_boundary_scratch_reset (&scratch);
        }
        return 0;
      }
    }
  }
  forest->set_from = forest_from;
  forest->t8code_data = data_temp;
  if (own_scratch) {
    balance_scratch = NULL;
    t8_forest_boundary_scratch_reset (&scratch);
  }
  return 1;
}

T8_EXTERN_C_END ();

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

/* The ghost layer of a given width (t8_forest_set_ghost_overlap) against a brute force reference: all leaves of all
 * processes are gathered, the elements of the overlap level are enumerated in a uniform forest on each process, and
 * their vertex neighborhoods are found by comparing vertex coordinates, independently of the neighbor computations
 * of the forest. The overlap ghosts of a process are the leaves of other processes whose ancestor of the overlap
 * level lies within the given number of vertex steps from an ancestor of a local leaf. Also checked: the overlap
 * contains the vertex ghost layer, and ghost data exchange delivers the data of the ghosts. */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <test/t8_gtest_macros.hxx>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

/** A leaf as (global tree, level, linear id at its level). */
using t8_test_leaf = std::tuple<t8_gloidx_t, int, t8_linearidx_t>;

/** Gather the leaves of all processes, with their ranks. */
static std::vector<std::pair<t8_test_leaf, int>>
t8_test_gather_leaves (const t8_forest_t forest)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  std::vector<int64_t> local;
  for (t8_locidx_t ltree = 0; ltree < t8_forest_get_num_local_trees (forest); ltree++) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltree);
    for (t8_locidx_t ileaf = 0; ileaf < t8_forest_get_tree_num_leaf_elements (forest, ltree); ileaf++) {
      const t8_element_t *leaf = t8_forest_get_leaf_element_in_tree (forest, ltree, ileaf);
      const int level = scheme->element_get_level (eclass, leaf);
      local.push_back (t8_forest_global_tree_id (forest, ltree));
      local.push_back (level);
      local.push_back ((int64_t) scheme->element_get_linear_id (eclass, leaf, level));
    }
  }
  const sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  int size, mpiret = sc_MPI_Comm_size (comm, &size);
  SC_CHECK_MPI (mpiret);
  int count = local.size ();
  std::vector<int> counts (size), displs (size + 1, 0);
  mpiret = sc_MPI_Allgather (&count, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  for (int rank = 0; rank < size; rank++) {
    displs[rank + 1] = displs[rank] + counts[rank];
  }
  std::vector<int64_t> all (displs[size]);
  mpiret = sc_MPI_Allgatherv (local.data (), count, sc_MPI_LONG_LONG_INT, all.data (), counts.data (), displs.data (),
                              sc_MPI_LONG_LONG_INT, comm);
  SC_CHECK_MPI (mpiret);
  std::vector<std::pair<t8_test_leaf, int>> leaves;
  for (int rank = 0; rank < size; rank++) {
    for (int i = displs[rank]; i < displs[rank + 1]; i += 3) {
      leaves.push_back ({ { all[i], (int) all[i + 1], (t8_linearidx_t) all[i + 2] }, rank });
    }
  }
  return leaves;
}

/** The ghosts of a forest as leaves. */
static std::set<t8_test_leaf>
t8_test_ghosts (const t8_forest_t forest)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  std::set<t8_test_leaf> ghosts;
  for (t8_locidx_t gtree = 0; gtree < t8_forest_ghost_num_trees (forest); gtree++) {
    const t8_eclass_t eclass = t8_forest_ghost_get_tree_class (forest, gtree);
    for (t8_locidx_t i = 0; i < t8_forest_ghost_tree_num_leaf_elements (forest, gtree); i++) {
      const t8_element_t *ghost = t8_forest_ghost_get_leaf_element (forest, gtree, i);
      const int level = scheme->element_get_level (eclass, ghost);
      ghosts.insert ({ t8_forest_ghost_get_global_treeid (forest, gtree), level,
                       scheme->element_get_linear_id (eclass, ghost, level) });
    }
  }
  return ghosts;
}

/** The overlap ghosts of this process by brute force: the leaves of other processes that intersect an element of the
 * level within the number of layers from an element of the level that intersects a local leaf. */
static std::set<t8_test_leaf>
t8_test_reference_overlap (const t8_forest_t forest, const t8_cmesh_t cmesh_self, const int level, const int layers)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const int rank = forest->mpirank;
  const std::vector<std::pair<t8_test_leaf, int>> leaves = t8_test_gather_leaves (forest);

  /* The elements of the level and their vertices, in a uniform forest on this process */
  t8_cmesh_ref (cmesh_self);
  scheme->ref ();
  t8_forest_t uniform = t8_forest_new_uniform (cmesh_self, scheme, level, 0, sc_MPI_COMM_SELF);
  std::vector<std::pair<t8_gloidx_t, t8_linearidx_t>> elements;
  std::map<std::array<long long, 3>, std::vector<size_t>> vertex_to_elements;
  for (t8_locidx_t ltree = 0; ltree < t8_forest_get_num_local_trees (uniform); ltree++) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (uniform, ltree);
    for (t8_locidx_t i = 0; i < t8_forest_get_tree_num_leaf_elements (uniform, ltree); i++) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (uniform, ltree, i);
      const size_t index = elements.size ();
      elements.push_back ({ t8_forest_global_tree_id (uniform, ltree),
                            scheme->element_get_linear_id (eclass, element, level) });
      for (int corner = 0; corner < scheme->element_get_num_corners (eclass, element); corner++) {
        double coords[3] = { 0, 0, 0 };
        t8_forest_element_coordinate (uniform, ltree, element, corner, coords);
        std::array<long long, 3> key;
        for (int d = 0; d < 3; d++) {
          key[d] = std::llround (coords[d] * 1e8);
        }
        vertex_to_elements[key].push_back (index);
      }
    }
  }
  std::map<std::pair<t8_gloidx_t, t8_linearidx_t>, size_t> element_index;
  for (size_t i = 0; i < elements.size (); i++) {
    element_index[elements[i]] = i;
  }
  std::vector<std::vector<size_t>> adjacent (elements.size ());
  for (const auto &[key, indices] : vertex_to_elements) {
    for (const size_t a : indices) {
      for (const size_t b : indices) {
        if (a != b) {
          adjacent[a].push_back (b);
        }
      }
    }
  }
  t8_forest_unref (&uniform);

  /* The elements of the level that a leaf intersects: its ancestor, or its descendants if it is coarser */
  const auto intersected = [&] (const t8_test_leaf &leaf) {
    const auto &[gtreeid, leaf_level, id] = leaf;
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh_self, gtreeid);
    t8_element_t *element, *desc;
    scheme->element_new (eclass, 1, &element);
    scheme->element_new (eclass, 1, &desc);
    scheme->element_set_linear_id (eclass, element, leaf_level, id);
    std::vector<size_t> result;
    if (leaf_level >= level) {
      result.push_back (element_index.at ({ gtreeid, scheme->element_get_linear_id (eclass, element, level) }));
    }
    else {
      scheme->element_get_first_descendant (eclass, element, desc, level);
      const t8_linearidx_t first = scheme->element_get_linear_id (eclass, desc, level);
      scheme->element_get_last_descendant (eclass, element, desc, level);
      const t8_linearidx_t last = scheme->element_get_linear_id (eclass, desc, level);
      for (t8_linearidx_t i = first; i <= last; i++) {
        result.push_back (element_index.at ({ gtreeid, i }));
      }
    }
    scheme->element_destroy (eclass, 1, &element);
    scheme->element_destroy (eclass, 1, &desc);
    return result;
  };
  /* The region: the elements intersecting local leaves, grown by the layers */
  std::vector<bool> in_region (elements.size (), false);
  std::vector<size_t> front;
  for (const auto &[leaf, owner] : leaves) {
    if (owner == rank) {
      for (const size_t a : intersected (leaf)) {
        if (!in_region[a]) {
          in_region[a] = true;
          front.push_back (a);
        }
      }
    }
  }
  for (int layer = 0; layer < layers; layer++) {
    std::vector<size_t> next;
    for (const size_t a : front) {
      for (const size_t b : adjacent[a]) {
        if (!in_region[b]) {
          in_region[b] = true;
          next.push_back (b);
        }
      }
    }
    front = next;
  }
  std::set<t8_test_leaf> reference;
  for (const auto &[leaf, owner] : leaves) {
    if (owner != rank) {
      for (const size_t a : intersected (leaf)) {
        if (in_region[a]) {
          reference.insert (leaf);
          break;
        }
      }
    }
  }
  return reference;
}

/* Refine the elements with (tree + linear id) divisible by 3 down to the level in the user data */
static int
t8_test_overlap_adapt (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                       const t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme_c *scheme,
                       [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                       t8_element_t *elements[])
{
  const int max_level = *(const int *) t8_forest_get_user_data (forest);
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest_from, which_tree);
  const t8_linearidx_t id = scheme->element_get_linear_id (tree_class, elements[0], level);
  return level < max_level && (gtreeid + id) % 3 == 0;
}

/** The mesh of a test: a hypercube of one element class, or the hybrid cube if eclass is T8_ECLASS_COUNT. */
static t8_cmesh_t
t8_test_overlap_cmesh (const t8_eclass_t eclass, sc_MPI_Comm comm)
{
  return eclass == T8_ECLASS_COUNT ? t8_cmesh_new_hypercube_hybrid (comm, 0, 0)
                                   : t8_cmesh_new_hypercube (eclass, comm, 0, 0, 0);
}

struct forest_ghost_overlap: public testing::TestWithParam<t8_eclass_t>
{
 protected:
  void
  SetUp () override
  {
    eclass = GetParam ();
    int size;
    const int mpiret = sc_MPI_Comm_size (sc_MPI_COMM_WORLD, &size);
    SC_CHECK_MPI (mpiret);
    if (size == 1) {
      GTEST_SKIP () << "No ghost layer on a single process";
    }
    scheme = t8_scheme_new_standalone ();
    cmesh = t8_test_overlap_cmesh (eclass, sc_MPI_COMM_WORLD);
    cmesh_self = t8_test_overlap_cmesh (eclass, sc_MPI_COMM_SELF);
    dim = eclass == T8_ECLASS_COUNT ? 3 : t8_eclass_to_dimension[eclass];
  }

  void
  TearDown () override
  {
    if (cmesh != NULL) {
      t8_cmesh_unref (&cmesh);
      t8_cmesh_unref (&cmesh_self);
      scheme->unref ();
    }
  }

  /** An adapted and partitioned forest, uniform up to \a min_level. */
  t8_forest_t
  adapted_forest (const int min_level)
  {
    t8_cmesh_ref (cmesh);
    scheme->ref ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, min_level, 0, sc_MPI_COMM_WORLD);
    max_level = min_level + 2;
    t8_forest_t forest_adapt;
    t8_forest_init (&forest_adapt);
    t8_forest_set_user_data (forest_adapt, &max_level);
    t8_forest_set_adapt (forest_adapt, forest, t8_test_overlap_adapt, 1);
    t8_forest_set_partition (forest_adapt, NULL, 0);
    t8_forest_commit (forest_adapt);
    return forest_adapt;
  }

  /** A copy of \a forest with the overlap ghost layer. */
  t8_forest_t
  with_overlap (const t8_forest_t forest, const int level, const int layers)
  {
    t8_forest_ref (forest);
    t8_forest_t copy;
    t8_forest_init (&copy);
    t8_forest_set_copy (copy, forest);
    t8_forest_set_ghost_overlap (copy, level, layers);
    t8_forest_commit (copy);
    return copy;
  }

  /** The overlap contains the vertex ghost layer, and ghost data exchange works. */
  void
  check_vertex_ghosts_and_exchange (const t8_forest_t forest, const t8_forest_t overlap)
  {
    t8_forest_ref (forest);
    t8_forest_t vertex;
    t8_forest_init (&vertex);
    t8_forest_set_copy (vertex, forest);
    t8_forest_set_ghost (vertex, 1, T8_GHOST_VERTICES);
    t8_forest_commit (vertex);
    const std::set<t8_test_leaf> vertex_ghosts = t8_test_ghosts (vertex), overlap_ghosts = t8_test_ghosts (overlap);
    for (const t8_test_leaf &leaf : vertex_ghosts) {
      EXPECT_TRUE (overlap_ghosts.count (leaf) > 0) << "a vertex ghost is not an overlap ghost";
    }
    t8_forest_unref (&vertex);

    /* Exchange a code of each leaf and compare it with the code of the ghost */
    const t8_scheme *s = t8_forest_get_scheme (overlap);
    const auto code = [&] (const t8_gloidx_t gtreeid, const t8_eclass_t tree_class, const t8_element_t *element) {
      const int level = s->element_get_level (tree_class, element);
      return (int64_t) gtreeid * 1000003 + level * 1009 + (int64_t) s->element_get_linear_id (tree_class, element, level);
    };
    const t8_locidx_t num_local = t8_forest_get_local_num_leaf_elements (overlap);
    const t8_locidx_t num_ghosts = t8_forest_get_num_ghosts (overlap);
    sc_array_t *data = sc_array_new_count (sizeof (int64_t), num_local + num_ghosts);
    t8_locidx_t index = 0;
    for (t8_locidx_t ltree = 0; ltree < t8_forest_get_num_local_trees (overlap); ltree++) {
      const t8_eclass_t tree_class = t8_forest_get_tree_class (overlap, ltree);
      for (t8_locidx_t i = 0; i < t8_forest_get_tree_num_leaf_elements (overlap, ltree); i++, index++) {
        *(int64_t *) sc_array_index (data, index) = code (t8_forest_global_tree_id (overlap, ltree), tree_class,
                                                         t8_forest_get_leaf_element_in_tree (overlap, ltree, i));
      }
    }
    /* The owners of the ghosts are known (deal.II looks them up for its ghost cells) */
    std::map<t8_test_leaf, int> owner_of;
    for (const auto &[leaf, owner] : t8_test_gather_leaves (overlap)) {
      owner_of[leaf] = owner;
    }
    for (t8_locidx_t gtree = 0; gtree < t8_forest_ghost_num_trees (overlap); gtree++) {
      const t8_eclass_t tree_class = t8_forest_ghost_get_tree_class (overlap, gtree);
      const t8_gloidx_t gtreeid = t8_forest_ghost_get_global_treeid (overlap, gtree);
      for (t8_locidx_t i = 0; i < t8_forest_ghost_tree_num_leaf_elements (overlap, gtree); i++) {
        t8_element_t *ghost = t8_forest_ghost_get_leaf_element (overlap, gtree, i);
        const int level = s->element_get_level (tree_class, ghost);
        EXPECT_EQ (t8_forest_element_find_owner (overlap, gtreeid, ghost, tree_class),
                   owner_of.at ({ gtreeid, level, s->element_get_linear_id (tree_class, ghost, level) }));
      }
    }

    t8_forest_ghost_exchange_data (overlap, data);
    for (t8_locidx_t gtree = 0; gtree < t8_forest_ghost_num_trees (overlap); gtree++) {
      const t8_eclass_t tree_class = t8_forest_ghost_get_tree_class (overlap, gtree);
      for (t8_locidx_t i = 0; i < t8_forest_ghost_tree_num_leaf_elements (overlap, gtree); i++, index++) {
        EXPECT_EQ (*(int64_t *) sc_array_index (data, index),
                   code (t8_forest_ghost_get_global_treeid (overlap, gtree), tree_class,
                         t8_forest_ghost_get_leaf_element (overlap, gtree, i)));
      }
    }
    EXPECT_EQ (index, num_local + num_ghosts);
    sc_array_destroy (data);
  }

  t8_eclass_t eclass;
  int dim;
  int max_level;
  const t8_scheme *scheme = NULL;
  t8_cmesh_t cmesh = NULL, cmesh_self = NULL;
};

TEST_P (forest_ghost_overlap, brute_force)
{
  const int min_level = dim == 3 ? 1 : 2;
  t8_forest_t forest = adapted_forest (min_level);
  for (int level = 0; level <= min_level; level++) {
    for (int layers = 1; layers <= 2; layers++) {
      t8_forest_t overlap = with_overlap (forest, level, layers);
      const std::set<t8_test_leaf> ghosts = t8_test_ghosts (overlap);
      const std::set<t8_test_leaf> reference = t8_test_reference_overlap (overlap, cmesh_self, level, layers);
      EXPECT_EQ (ghosts, reference) << "level " << level << " layers " << layers << ": " << ghosts.size ()
                                    << " ghosts, " << reference.size () << " expected";
      check_vertex_ghosts_and_exchange (forest, overlap);
      t8_forest_unref (&overlap);
    }
  }
  t8_forest_unref (&forest);
}

TEST_P (forest_ghost_overlap, coarse_leaves_and_empty_processes)
{
  /* Leaves coarser than the overlap level, and at level 0 processes without leaves */
  for (const int min_level : { 0, 1 }) {
    t8_forest_t forest = adapted_forest (min_level);
    for (const int layers : { 1, 2 }) {
      const int level = min_level + 2;
      t8_forest_t overlap = with_overlap (forest, level, layers);
      const std::set<t8_test_leaf> ghosts = t8_test_ghosts (overlap);
      const std::set<t8_test_leaf> reference = t8_test_reference_overlap (overlap, cmesh_self, level, layers);
      EXPECT_EQ (ghosts, reference) << "level " << level << " layers " << layers << ": " << ghosts.size ()
                                    << " ghosts, " << reference.size () << " expected";
      check_vertex_ghosts_and_exchange (forest, overlap);
      t8_forest_unref (&overlap);
    }
    t8_forest_unref (&forest);
  }
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_overlap, forest_ghost_overlap,
                          ::testing::Values (T8_ECLASS_QUAD, T8_ECLASS_TRIANGLE, T8_ECLASS_HEX, T8_ECLASS_TET,
                                             T8_ECLASS_PRISM, T8_ECLASS_PYRAMID, T8_ECLASS_COUNT),
                          [] (const testing::TestParamInfo<t8_eclass_t> &info) {
                            return info.param == T8_ECLASS_COUNT ? std::string ("hybrid")
                                                                 : std::string (t8_eclass_to_string[info.param]);
                          });

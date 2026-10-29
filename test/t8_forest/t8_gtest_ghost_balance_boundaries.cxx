/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_gtest_ghost_balance_boundaries.cxx
 * Tests face, edge and vertex ghosts and balance of unbalanced forests of the standalone scheme on a single tree.
 * The reference is a brute force computation: all leaves are gathered on all ranks, and two leaves are neighbors
 * if they touch with a boundary of at least the dimension of the ghost/balance type. Since the leaves stem from
 * refining a conforming mesh, the intersection of two leaves is a boundary of the finer one, spanned by its corners
 * that lie in the coarser one.
 */

#include <gtest/gtest.h>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_schemes/t8_standalone/t8_standalone.hxx>
#include <t8_schemes/t8_standalone/t8_standalone_elements.hxx>
#include <test/t8_gtest_macros.hxx>
#include <array>
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

/* A leaf in a form independent of the eclass */
struct t8_test_leaf
{
  std::array<int64_t, 3> coords = { 0, 0, 0 };
  int level = 0;
  int type = 0;
  bool
  operator< (const t8_test_leaf &other) const
  {
    return std::tie (level, type, coords) < std::tie (other.level, other.type, other.coords);
  }
};

/* The geometry of the types of an eclass: element x lies in the type t element at anchor a with length h iff
 * 0 <= x_i - a_i <= h and for every equation e = (i, j): x_i - a_i >= x_j - a_j if bit e of t is 0, else <=. */
struct t8_test_shape
{
  int dim = 0;
  int maxlevel = 0;
  std::vector<std::array<int, 2>> equations;
  std::vector<std::vector<int>> corners; /* per type the corners as bitsets of the cube */

  int64_t
  length (const t8_test_leaf &leaf) const
  {
    return (int64_t) 1 << (maxlevel - leaf.level);
  }

  bool
  contains (const t8_test_leaf &leaf, const std::array<int64_t, 3> &point) const
  {
    const int64_t h = length (leaf);
    std::array<int64_t, 3> rel;
    for (int idim = 0; idim < dim; idim++) {
      rel[idim] = point[idim] - leaf.coords[idim];
      if (rel[idim] < 0 || rel[idim] > h) {
        return false;
      }
    }
    for (size_t ieq = 0; ieq < equations.size (); ieq++) {
      const int i = equations[ieq][0], j = equations[ieq][1];
      if (((leaf.type >> ieq) & 1) ? rel[i] > rel[j] : rel[i] < rel[j]) {
        return false;
      }
    }
    return true;
  }

  /* The dimension of the intersection of two leaves, -1 if they do not touch */
  int
  touch_dim (const t8_test_leaf &a, const t8_test_leaf &b) const
  {
    const t8_test_leaf &fine = a.level >= b.level ? a : b;
    const t8_test_leaf &coarse = a.level >= b.level ? b : a;
    const int64_t h = length (fine);
    /* Corners of the fine leaf in the coarse leaf, relative to the first one in units of h */
    std::vector<std::array<double, 3>> vectors;
    std::array<int64_t, 3> first;
    bool found = false;
    for (const int corner : corners[fine.type]) {
      std::array<int64_t, 3> point = fine.coords;
      for (int idim = 0; idim < dim; idim++) {
        point[idim] += ((corner >> idim) & 1) * h;
      }
      if (!contains (coarse, point)) {
        continue;
      }
      if (!found) {
        first = point;
        found = true;
      }
      else {
        std::array<double, 3> v = { 0, 0, 0 };
        for (int idim = 0; idim < dim; idim++) {
          v[idim] = (double) ((point[idim] - first[idim]) / h);
        }
        vectors.push_back (v);
      }
    }
    if (!found) {
      return -1;
    }
    /* Rank of the difference vectors by Gaussian elimination, the entries are small integers */
    int rank = 0;
    for (int col = 0; col < dim && rank < (int) vectors.size (); col++) {
      size_t pivot = rank;
      while (pivot < vectors.size () && std::fabs (vectors[pivot][col]) < 0.5) {
        pivot++;
      }
      if (pivot == vectors.size ()) {
        continue;
      }
      std::swap (vectors[rank], vectors[pivot]);
      for (size_t row = 0; row < vectors.size (); row++) {
        if ((int) row != rank && std::fabs (vectors[row][col]) > 1e-9) {
          const double factor = vectors[row][col] / vectors[rank][col];
          for (int c = 0; c < dim; c++) {
            vectors[row][c] -= factor * vectors[rank][c];
          }
        }
      }
      rank++;
    }
    return rank;
  }
};

template <t8_eclass TEclass>
static t8_test_shape
t8_test_make_shape ()
{
  t8_test_shape shape;
  shape.dim = T8_ELEMENT_DIM[TEclass];
  shape.maxlevel = T8_ELEMENT_MAXLEVEL[TEclass];
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
    for (int ieq = 0; ieq < T8_ELEMENT_NUM_EQUATIONS[TEclass]; ieq++) {
      shape.equations.push_back ({ t8_type_edge_equations<TEclass>[ieq][0], t8_type_edge_equations<TEclass>[ieq][1] });
    }
    for (int type = 0; type < (1 << T8_ELEMENT_NUM_EQUATIONS[TEclass]); type++) {
      std::vector<int> corners;
      for (int icorner = 0; icorner < T8_ELEMENT_NUM_CORNERS[TEclass]; icorner++) {
        int corner = 0;
        bool valid = true;
        for (int idim = 0; idim < shape.dim; idim++) {
          const int bit = t8_type_vertex_dim_to_binary<TEclass>[type][icorner][idim];
          valid = valid && bit != -1;
          corner |= (bit == 1) << idim;
        }
        if (valid) {
          corners.push_back (corner);
        }
      }
      shape.corners.push_back (corners);
    }
  }
  else {
    std::vector<int> corners;
    for (int corner = 0; corner < (1 << shape.dim); corner++) {
      corners.push_back (corner);
    }
    shape.corners.push_back (corners);
  }
  return shape;
}

template <t8_eclass TEclass>
static t8_test_leaf
t8_test_leaf_from_element_impl (const t8_element_t *element)
{
  const t8_standalone_element<TEclass> *el = (const t8_standalone_element<TEclass> *) element;
  t8_test_leaf leaf;
  for (int idim = 0; idim < T8_ELEMENT_DIM[TEclass]; idim++) {
    leaf.coords[idim] = el->coords[idim];
  }
  leaf.level = el->level;
  if constexpr (T8_ELEMENT_NUM_EQUATIONS[TEclass]) {
    leaf.type = el->type;
  }
  return leaf;
}

#define T8_TEST_ECLASS_DISPATCH(eclass, function, ...) \
  switch (eclass) { \
  case T8_ECLASS_LINE: \
    return function<T8_ECLASS_LINE> (__VA_ARGS__); \
  case T8_ECLASS_QUAD: \
    return function<T8_ECLASS_QUAD> (__VA_ARGS__); \
  case T8_ECLASS_TRIANGLE: \
    return function<T8_ECLASS_TRIANGLE> (__VA_ARGS__); \
  case T8_ECLASS_HEX: \
    return function<T8_ECLASS_HEX> (__VA_ARGS__); \
  case T8_ECLASS_TET: \
    return function<T8_ECLASS_TET> (__VA_ARGS__); \
  case T8_ECLASS_PRISM: \
    return function<T8_ECLASS_PRISM> (__VA_ARGS__); \
  case T8_ECLASS_PYRAMID: \
    return function<T8_ECLASS_PYRAMID> (__VA_ARGS__); \
  default: \
    SC_ABORT_NOT_REACHED (); \
  }

static t8_test_shape
t8_test_make_shape (const t8_eclass_t eclass)
{
  T8_TEST_ECLASS_DISPATCH (eclass, t8_test_make_shape);
}

static t8_test_leaf
t8_test_leaf_from_element (const t8_eclass_t eclass, const t8_element_t *element)
{
  T8_TEST_ECLASS_DISPATCH (eclass, t8_test_leaf_from_element_impl, element);
}

/* The smallest dimension of a touching boundary that makes two leaves neighbors for a ghost/balance type */
static int
t8_test_min_dim (const t8_ghost_type_t type, const int dim)
{
  switch (type) {
  case T8_GHOST_FACES:
    return SC_MAX (dim - 1, 0);
  case T8_GHOST_EDGES:
    return SC_MAX (dim - 2, 0);
  default:
    return 0;
  }
}

/* The local leaves of a single tree forest */
static std::vector<t8_test_leaf>
t8_test_local_leaves (t8_forest_t forest)
{
  std::vector<t8_test_leaf> leaves;
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ielem++) {
      leaves.push_back (t8_test_leaf_from_element (eclass, t8_forest_get_leaf_element_in_tree (forest, itree, ielem)));
    }
  }
  return leaves;
}

/* All leaves of the forest with their owners, gathered on all ranks */
static void
t8_test_all_leaves (t8_forest_t forest, std::vector<t8_test_leaf> &leaves, std::vector<int> &owners)
{
  const sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  int mpisize;
  int mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  /* Pack each leaf as 5 integers */
  std::vector<int64_t> send;
  for (const t8_test_leaf &leaf : t8_test_local_leaves (forest)) {
    send.insert (send.end (), { leaf.coords[0], leaf.coords[1], leaf.coords[2], (int64_t) leaf.level,
                                (int64_t) leaf.type });
  }
  int send_count = send.size ();
  std::vector<int> counts (mpisize), offsets (mpisize, 0);
  mpiret = sc_MPI_Allgather (&send_count, 1, sc_MPI_INT, counts.data (), 1, sc_MPI_INT, comm);
  SC_CHECK_MPI (mpiret);
  for (int rank = 1; rank < mpisize; rank++) {
    offsets[rank] = offsets[rank - 1] + counts[rank - 1];
  }
  std::vector<int64_t> recv (offsets[mpisize - 1] + counts[mpisize - 1]);
  mpiret = sc_MPI_Allgatherv (send.data (), send_count, sc_MPI_LONG_LONG_INT, recv.data (), counts.data (),
                              offsets.data (), sc_MPI_LONG_LONG_INT, comm);
  SC_CHECK_MPI (mpiret);
  leaves.clear ();
  owners.clear ();
  for (int rank = 0; rank < mpisize; rank++) {
    for (int ientry = offsets[rank]; ientry < offsets[rank] + counts[rank]; ientry += 5) {
      t8_test_leaf leaf;
      leaf.coords = { recv[ientry], recv[ientry + 1], recv[ientry + 2] };
      leaf.level = recv[ientry + 3];
      leaf.type = recv[ientry + 4];
      leaves.push_back (leaf);
      owners.push_back (rank);
    }
  }
}

/* Adapt callback: refine all elements that contain a point up to a deep level, and pseudo randomly some other
 * elements up to a coarse level, which gives an unbalanced forest. */
struct t8_test_adapt_data
{
  t8_test_shape shape;
  std::array<int64_t, 3> point;
  int deep_level;
  int random_level;
};

static int
t8_test_adapt_unbalanced (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from,
                          [[maybe_unused]] t8_locidx_t which_tree, t8_eclass_t tree_class,
                          [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                          [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                          t8_element_t *elements[])
{
  const t8_test_adapt_data *data = (const t8_test_adapt_data *) t8_forest_get_user_data (forest);
  const t8_test_leaf leaf = t8_test_leaf_from_element (tree_class, elements[0]);
  if (leaf.level < data->deep_level && data->shape.contains (leaf, data->point)) {
    return 1;
  }
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_linearidx_t id = scheme->element_get_linear_id (tree_class, elements[0], level);
  return level < data->random_level && (id * 7919 + level) % 3 == 0;
}

struct forest_ghost_balance_boundaries: public testing::TestWithParam<std::tuple<t8_eclass_t, t8_ghost_type_t>>
{
 protected:
  void
  SetUp () override
  {
    eclass = std::get<0> (GetParam ());
    type = std::get<1> (GetParam ());
    data.shape = t8_test_make_shape (eclass);
    const int dim = data.shape.dim;
    /* A point in the interior of the root with odd coordinates at level deep_level */
#if T8_TEST_LEVEL_INT == 0
    data.deep_level = dim == 3 ? 4 : 6;
    data.random_level = 2;
#elif T8_TEST_LEVEL_INT == 1
    data.deep_level = dim == 3 ? 3 : 5;
    data.random_level = 2;
#else
    data.deep_level = dim == 3 ? 3 : 4;
    data.random_level = 1;
#endif
    const int64_t h = (int64_t) 1 << (data.shape.maxlevel - data.deep_level);
    const int64_t root_len = (int64_t) 1 << data.shape.maxlevel;
    data.point = { 0, 0, 0 };
    for (int idim = 0; idim < dim; idim++) {
      /* Decreasing coordinates lie in the root of type 0 of all simplices, prisms and pyramids */
      data.point[idim] = (root_len / 2 + root_len / 8 - idim * (root_len / 8)) / h * h + h * (idim == 0);
    }
    t8_cmesh_t cmesh = t8_cmesh_new_from_class (eclass, sc_MPI_COMM_WORLD);
    const t8_scheme *scheme = t8_scheme_new_standalone ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 1, 0, sc_MPI_COMM_WORLD);
    forest = t8_forest_new_adapt (forest, t8_test_adapt_unbalanced, 1, 0, &data);
    /* Partition and create the ghost layer of the type */
    t8_forest_init (&unbalanced);
    t8_forest_set_partition (unbalanced, forest, 0);
    t8_forest_set_ghost (unbalanced, 1, type);
    t8_forest_commit (unbalanced);
  }

  void
  TearDown () override
  {
    t8_forest_unref (&unbalanced);
  }

  /* Whether each leaf of the forest has only neighbors of the type whose level differs by at most one
   * (logical and over all ranks) */
  int
  is_balanced_brute_force (t8_forest_t forest)
  {
    std::vector<t8_test_leaf> all_leaves;
    std::vector<int> owners;
    t8_test_all_leaves (forest, all_leaves, owners);
    const int min_dim = t8_test_min_dim (type, data.shape.dim);
    int balanced = 1;
    for (const t8_test_leaf &leaf : t8_test_local_leaves (forest)) {
      for (const t8_test_leaf &other : all_leaves) {
        if (std::abs (leaf.level - other.level) > 1 && data.shape.touch_dim (leaf, other) >= min_dim) {
          balanced = 0;
        }
      }
    }
    int balanced_global;
    const int mpiret
      = sc_MPI_Allreduce (&balanced, &balanced_global, 1, sc_MPI_INT, sc_MPI_LAND, t8_forest_get_mpicomm (forest));
    SC_CHECK_MPI (mpiret);
    return balanced_global;
  }

  t8_eclass_t eclass;
  t8_ghost_type_t type;
  t8_test_adapt_data data;
  t8_forest_t unbalanced;
};

/* The ghost layer consists of exactly the remote leaves that touch a local leaf with a boundary of the type */
TEST_P (forest_ghost_balance_boundaries, ghost_layer)
{
  std::vector<t8_test_leaf> all_leaves;
  std::vector<int> owners;
  t8_test_all_leaves (unbalanced, all_leaves, owners);
  const int min_dim = t8_test_min_dim (type, data.shape.dim);
  int mpirank;
  const int mpiret = sc_MPI_Comm_rank (t8_forest_get_mpicomm (unbalanced), &mpirank);
  SC_CHECK_MPI (mpiret);

  std::set<t8_test_leaf> expected;
  const std::vector<t8_test_leaf> local_leaves = t8_test_local_leaves (unbalanced);
  for (size_t iother = 0; iother < all_leaves.size (); iother++) {
    if (owners[iother] == mpirank) {
      continue;
    }
    for (const t8_test_leaf &leaf : local_leaves) {
      if (data.shape.touch_dim (leaf, all_leaves[iother]) >= min_dim) {
        expected.insert (all_leaves[iother]);
        break;
      }
    }
  }

  std::set<t8_test_leaf> ghosts;
  for (t8_locidx_t itree = 0; itree < t8_forest_ghost_num_trees (unbalanced); itree++) {
    for (t8_locidx_t ielem = 0; ielem < t8_forest_ghost_tree_num_leaf_elements (unbalanced, itree); ielem++) {
      ghosts.insert (t8_test_leaf_from_element (eclass, t8_forest_ghost_get_leaf_element (unbalanced, itree, ielem)));
    }
  }
  EXPECT_EQ (ghosts.size (), expected.size ());
  for (const t8_test_leaf &leaf : expected) {
    EXPECT_TRUE (ghosts.count (leaf)) << "Missing ghost at level " << leaf.level << " type " << leaf.type
                                      << " anchor " << leaf.coords[0] << " " << leaf.coords[1] << " "
                                      << leaf.coords[2];
  }
  for (const t8_test_leaf &leaf : ghosts) {
    EXPECT_TRUE (expected.count (leaf)) << "Ghost that does not touch a local leaf at level " << leaf.level;
  }
}

/* Balance makes the forest balanced for the type, and t8_forest_is_balanced agrees */
TEST_P (forest_ghost_balance_boundaries, balance)
{
  /* The test forest is not balanced */
  EXPECT_FALSE (is_balanced_brute_force (unbalanced));

  t8_forest_ref (unbalanced);
  t8_forest_t balanced;
  t8_forest_init (&balanced);
  t8_forest_set_balance (balanced, unbalanced, 0);
  t8_forest_set_balance_type (balanced, type);
  t8_forest_set_ghost (balanced, 1, type);
  t8_forest_commit (balanced);

  EXPECT_TRUE (is_balanced_brute_force (balanced));
  EXPECT_TRUE (t8_forest_is_balanced (balanced));

  if (type != T8_GHOST_FACES) {
    /* Balancing more neighbors refines at least as much as face balance */
    t8_forest_ref (unbalanced);
    t8_forest_t face_balanced;
    t8_forest_init (&face_balanced);
    t8_forest_set_balance (face_balanced, unbalanced, 0);
    t8_forest_commit (face_balanced);
    EXPECT_GE (t8_forest_get_global_num_leaf_elements (balanced),
               t8_forest_get_global_num_leaf_elements (face_balanced));
    t8_forest_unref (&face_balanced);
  }

  /* Balancing again does not change the forest */
  t8_forest_ref (balanced);
  t8_forest_t balanced_again;
  t8_forest_init (&balanced_again);
  t8_forest_set_balance (balanced_again, balanced, 0);
  t8_forest_set_balance_type (balanced_again, type);
  t8_forest_commit (balanced_again);
  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (balanced_again),
             t8_forest_get_global_num_leaf_elements (balanced));
  t8_forest_unref (&balanced_again);
  t8_forest_unref (&balanced);
}

auto print_eclass_ghost_type = [] (const testing::TestParamInfo<std::tuple<t8_eclass_t, t8_ghost_type_t>> &info) {
  static const char *type_names[] = { "none", "faces", "edges", "vertices" };
  return std::string (t8_eclass_to_string[std::get<0> (info.param)]) + "_" + type_names[std::get<1> (info.param)];
};

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_balance_boundaries, forest_ghost_balance_boundaries,
                          ::testing::Combine (::testing::Range (T8_ECLASS_LINE, T8_ECLASS_COUNT),
                                              ::testing::Values (T8_GHOST_FACES, T8_GHOST_EDGES, T8_GHOST_VERTICES)),
                          print_eclass_ghost_type);

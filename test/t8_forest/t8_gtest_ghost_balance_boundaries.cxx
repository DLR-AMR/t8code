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
 * Tests face, edge and vertex ghosts and balance of unbalanced forests of the standalone scheme on single trees and
 * on meshes of several trees. The reference is a brute force computation: all leaves are gathered on all ranks, and
 * two leaves are neighbors if they touch with a boundary of at least the dimension of the ghost/balance type. Since
 * the leaves stem from refining a conforming mesh, the intersection of two leaves is a boundary of the finer one,
 * spanned by its corners that lie in the coarser one. The trees of the meshes are affine images of their reference
 * elements, so the corners of a leaf are mapped into the tree of another leaf by the affine maps of the trees.
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
#include <functional>
#include <map>
#include <cmath>
#include <set>
#include <tuple>
#include <vector>

/* A leaf in a form independent of the eclass */
struct t8_test_leaf
{
  t8_gloidx_t gtreeid = 0;
  std::array<int64_t, 3> coords = { 0, 0, 0 };
  int level = 0;
  int type = 0;
  bool
  operator< (const t8_test_leaf &other) const
  {
    return std::tie (gtreeid, level, type, coords) < std::tie (other.gtreeid, other.level, other.type, other.coords);
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

  int64_t
  root_length () const
  {
    return (int64_t) 1 << maxlevel;
  }

  template <typename TCoord>
  bool
  contains (const t8_test_leaf &leaf, const std::array<TCoord, 3> &point) const
  {
    const int64_t h = length (leaf);
    std::array<TCoord, 3> rel;
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

  /* The corners of a leaf in the coordinates of its tree */
  std::vector<std::array<int64_t, 3>>
  leaf_corners (const t8_test_leaf &leaf) const
  {
    std::vector<std::array<int64_t, 3>> points;
    for (const int corner : corners[leaf.type]) {
      std::array<int64_t, 3> point = leaf.coords;
      for (int idim = 0; idim < dim; idim++) {
        point[idim] += ((corner >> idim) & 1) * length (leaf);
      }
      points.push_back (point);
    }
    return points;
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
t8_test_leaf_from_element (const t8_eclass_t eclass, const t8_gloidx_t gtreeid, const t8_element_t *element)
{
  auto leaf_from_element = [&] () { T8_TEST_ECLASS_DISPATCH (eclass, t8_test_leaf_from_element_impl, element); };
  t8_test_leaf leaf = leaf_from_element ();
  leaf.gtreeid = gtreeid;
  return leaf;
}

/* The trees of a mesh with their affine maps x -> A x / L + b from the coordinates of the tree (with root length L)
 * to the first dim physical coordinates. */
struct t8_test_mesh
{
  int dim = 0;
  std::vector<t8_eclass_t> eclasses;
  std::map<t8_eclass_t, t8_test_shape> shapes;
  std::vector<std::array<std::array<double, 3>, 3>> matrices, inverses;
  std::vector<std::array<double, 3>> offsets;

  const t8_test_shape &
  shape (const t8_test_leaf &leaf) const
  {
    return shapes.at (eclasses[leaf.gtreeid]);
  }

  std::array<double, 3>
  to_physical (const t8_gloidx_t gtreeid, const std::array<int64_t, 3> &point) const
  {
    const double root_length = shapes.at (eclasses[gtreeid]).root_length ();
    std::array<double, 3> physical = offsets[gtreeid];
    for (int idim = 0; idim < dim; idim++) {
      for (int jdim = 0; jdim < dim; jdim++) {
        physical[idim] += matrices[gtreeid][idim][jdim] * point[jdim] / root_length;
      }
    }
    return physical;
  }

  /* The coordinates of a physical point in a tree, and whether they are integers */
  bool
  to_tree (const t8_gloidx_t gtreeid, const std::array<double, 3> &physical, std::array<int64_t, 3> &point) const
  {
    const double root_length = shapes.at (eclasses[gtreeid]).root_length ();
    point = { 0, 0, 0 };
    for (int idim = 0; idim < dim; idim++) {
      double coord = 0;
      for (int jdim = 0; jdim < dim; jdim++) {
        coord += inverses[gtreeid][idim][jdim] * (physical[jdim] - offsets[gtreeid][jdim]) * root_length;
      }
      point[idim] = std::llround (coord);
      if (std::fabs (coord - point[idim]) > 1e-3) {
        return false;
      }
    }
    return true;
  }

  /* Whether a leaf contains a physical point, which need not be a corner of the refinement */
  bool
  contains (const t8_test_leaf &leaf, const std::array<double, 3> &physical) const
  {
    const double root_length = shape (leaf).root_length ();
    std::array<double, 3> point = { 0, 0, 0 };
    for (int idim = 0; idim < dim; idim++) {
      for (int jdim = 0; jdim < dim; jdim++) {
        point[idim] += inverses[leaf.gtreeid][idim][jdim] * (physical[jdim] - offsets[leaf.gtreeid][jdim]) * root_length;
      }
    }
    return shape (leaf).contains (leaf, point);
  }

  /* The dimension of the intersection of two leaves, -1 if they do not touch */
  int
  touch_dim (const t8_test_leaf &a, const t8_test_leaf &b) const
  {
    const t8_test_leaf &fine = a.level >= b.level ? a : b;
    const t8_test_leaf &coarse = a.level >= b.level ? b : a;
    /* Corners of the fine leaf in the coarse leaf, relative to the first one */
    std::vector<std::array<double, 3>> vectors;
    std::array<int64_t, 3> first;
    bool found = false;
    for (const auto &corner : shape (fine).leaf_corners (fine)) {
      std::array<int64_t, 3> point;
      if (!to_tree (coarse.gtreeid, to_physical (fine.gtreeid, corner), point)
          || !shape (coarse).contains (coarse, point)) {
        continue;
      }
      if (!found) {
        first = point;
        found = true;
      }
      else {
        std::array<double, 3> v = { 0, 0, 0 };
        for (int idim = 0; idim < dim; idim++) {
          v[idim] = (double) (point[idim] - first[idim]) / ((int64_t) 1 << (shape (coarse).maxlevel - fine.level));
        }
        vectors.push_back (v);
      }
    }
    if (!found) {
      return -1;
    }
    /* Rank of the difference vectors by Gaussian elimination */
    int rank = 0;
    for (int col = 0; col < dim && rank < (int) vectors.size (); col++) {
      size_t pivot = rank;
      while (pivot < vectors.size () && std::fabs (vectors[pivot][col]) < 1e-9) {
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

/* Invert a dim x dim matrix by Gauss-Jordan elimination */
static std::array<std::array<double, 3>, 3>
t8_test_invert (std::array<std::array<double, 3>, 3> matrix, const int dim)
{
  std::array<std::array<double, 3>, 3> inverse = {};
  for (int idim = 0; idim < dim; idim++) {
    inverse[idim][idim] = 1;
  }
  for (int col = 0; col < dim; col++) {
    int pivot = col;
    for (int row = col + 1; row < dim; row++) {
      if (std::fabs (matrix[row][col]) > std::fabs (matrix[pivot][col])) {
        pivot = row;
      }
    }
    std::swap (matrix[col], matrix[pivot]);
    std::swap (inverse[col], inverse[pivot]);
    const double diagonal = matrix[col][col];
    for (int c = 0; c < dim; c++) {
      matrix[col][c] /= diagonal;
      inverse[col][c] /= diagonal;
    }
    for (int row = 0; row < dim; row++) {
      if (row != col) {
        const double factor = matrix[row][col];
        for (int c = 0; c < dim; c++) {
          matrix[row][c] -= factor * matrix[col][c];
          inverse[row][c] -= factor * inverse[col][c];
        }
      }
    }
  }
  return inverse;
}

/* The affine maps of the trees of a replicated cmesh, from the tree vertices */
static t8_test_mesh
t8_test_make_mesh (const t8_cmesh_t cmesh)
{
  t8_test_mesh mesh;
  for (t8_gloidx_t itree = 0; itree < t8_cmesh_get_num_trees (cmesh); itree++) {
    const t8_eclass_t eclass = t8_cmesh_get_tree_class (cmesh, itree);
    if (!mesh.shapes.count (eclass)) {
      mesh.shapes[eclass] = t8_test_make_shape (eclass);
    }
    const t8_test_shape &shape = mesh.shapes[eclass];
    mesh.dim = shape.dim;
    mesh.eclasses.push_back (eclass);
    const double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
    /* The root corners are the tree vertices. Corner 0 is the origin, and we take dim corners whose cube vertices
     * are a chain 1, 3, 7 (simplices and pyramids) or the unit vectors (cubes and prisms). */
    const std::vector<int> &corners = shape.corners[0];
    std::array<double, 3> offset = { vertices[0], vertices[1], vertices[2] };
    std::array<std::array<double, 3>, 3> unit_images = {}, units = {};
    int num_basis = 0;
    for (size_t icorner = 1; icorner < corners.size () && num_basis < shape.dim; icorner++) {
      /* Use the corner if its cube vertex is independent of the ones used so far */
      std::array<double, 3> unit = { 0, 0, 0 };
      for (int idim = 0; idim < shape.dim; idim++) {
        unit[idim] = (corners[icorner] >> idim) & 1;
      }
      std::array<std::array<double, 3>, 3> test = units;
      test[num_basis] = unit;
      bool independent = true;
      {
        /* Check the rank with a determinant-free elimination on a copy */
        std::array<std::array<double, 3>, 3> m = test;
        int rank = 0;
        for (int col = 0; col < shape.dim && rank <= num_basis; col++) {
          int pivot = rank;
          while (pivot <= num_basis && std::fabs (m[pivot][col]) < 1e-9) {
            pivot++;
          }
          if (pivot > num_basis) {
            continue;
          }
          std::swap (m[rank], m[pivot]);
          for (int row = 0; row <= num_basis; row++) {
            if (row != rank) {
              const double factor = m[row][col] / m[rank][col];
              for (int c = 0; c < shape.dim; c++) {
                m[row][c] -= factor * m[rank][c];
              }
            }
          }
          rank++;
        }
        independent = rank == num_basis + 1;
      }
      if (independent) {
        units[num_basis] = unit;
        for (int idim = 0; idim < shape.dim; idim++) {
          unit_images[num_basis][idim] = vertices[3 * icorner + idim] - offset[idim];
        }
        num_basis++;
      }
    }
    EXPECT_EQ (num_basis, shape.dim);
    /* A U = P with the columns u_i of U and p_i of P, so A = P U^-1 */
    const auto units_inverse = t8_test_invert (units, shape.dim);
    std::array<std::array<double, 3>, 3> matrix = {};
    for (int idim = 0; idim < shape.dim; idim++) {
      for (int jdim = 0; jdim < shape.dim; jdim++) {
        for (int k = 0; k < shape.dim; k++) {
          /* P[idim][k] = unit_images[k][idim], U^-1[k][jdim] = units_inverse transposed, since units are rows */
          matrix[idim][jdim] += unit_images[k][idim] * units_inverse[jdim][k];
        }
      }
    }
    mesh.matrices.push_back (matrix);
    mesh.inverses.push_back (t8_test_invert (matrix, shape.dim));
    mesh.offsets.push_back (offset);
    /* The tree is an affine image of its reference element */
    for (size_t icorner = 0; icorner < corners.size (); icorner++) {
      std::array<int64_t, 3> point = { 0, 0, 0 };
      for (int idim = 0; idim < shape.dim; idim++) {
        point[idim] = ((corners[icorner] >> idim) & 1) * shape.root_length ();
      }
      const auto physical = mesh.to_physical (itree, point);
      for (int idim = 0; idim < shape.dim; idim++) {
        EXPECT_NEAR (physical[idim], vertices[3 * icorner + idim], 1e-9) << "tree " << itree << " corner " << icorner;
      }
    }
  }
  return mesh;
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

/* The local leaves of a forest */
static std::vector<t8_test_leaf>
t8_test_local_leaves (t8_forest_t forest)
{
  std::vector<t8_test_leaf> leaves;
  for (t8_locidx_t itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, itree);
    const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_get_tree_num_leaf_elements (forest, itree); ielem++) {
      leaves.push_back (
        t8_test_leaf_from_element (eclass, gtreeid, t8_forest_get_leaf_element_in_tree (forest, itree, ielem)));
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
  /* Pack each leaf as 6 integers */
  std::vector<int64_t> send;
  for (const t8_test_leaf &leaf : t8_test_local_leaves (forest)) {
    send.insert (send.end (), { (int64_t) leaf.gtreeid, leaf.coords[0], leaf.coords[1], leaf.coords[2],
                                (int64_t) leaf.level, (int64_t) leaf.type });
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
    for (int ientry = offsets[rank]; ientry < offsets[rank] + counts[rank]; ientry += 6) {
      t8_test_leaf leaf;
      leaf.gtreeid = recv[ientry];
      leaf.coords = { recv[ientry + 1], recv[ientry + 2], recv[ientry + 3] };
      leaf.level = recv[ientry + 4];
      leaf.type = recv[ientry + 5];
      leaves.push_back (leaf);
      owners.push_back (rank);
    }
  }
}

/* Adapt callback: refine all elements that contain a point up to a deep level, and pseudo randomly some other
 * elements up to a coarse level, which gives an unbalanced forest. */
struct t8_test_adapt_data
{
  t8_test_mesh mesh;
  std::array<double, 3> point;
  int deep_level;
  int random_level;
};

static int
t8_test_adapt_unbalanced (t8_forest_t forest, [[maybe_unused]] t8_forest_t forest_from, t8_locidx_t which_tree,
                          t8_eclass_t tree_class, [[maybe_unused]] t8_locidx_t lelement_id, const t8_scheme *scheme,
                          [[maybe_unused]] const int is_family, [[maybe_unused]] const int num_elements,
                          t8_element_t *elements[])
{
  const t8_test_adapt_data *data = (const t8_test_adapt_data *) t8_forest_get_user_data (forest);
  const t8_test_leaf leaf
    = t8_test_leaf_from_element (tree_class, t8_forest_global_tree_id (forest_from, which_tree), elements[0]);
  if (leaf.level < data->deep_level && data->mesh.contains (leaf, data->point)) {
    return 1;
  }
  const int level = scheme->element_get_level (tree_class, elements[0]);
  const t8_linearidx_t id = scheme->element_get_linear_id (tree_class, elements[0], level);
  return level < data->random_level && (id * 7919 + level + leaf.gtreeid) % 3 == 0;
}

/* A test mesh: its name and a function that creates the replicated cmesh */
using t8_test_cmesh = std::pair<std::string, std::function<t8_cmesh_t ()>>;

struct forest_ghost_balance_boundaries: public testing::TestWithParam<std::tuple<t8_test_cmesh, t8_ghost_type_t>>
{
 protected:
  void
  SetUp () override
  {
    type = std::get<1> (GetParam ());
    t8_cmesh_t cmesh = std::get<0> (GetParam ()).second ();
    data.mesh = t8_test_make_mesh (cmesh);
    const int dim = data.mesh.dim;
    const int num_trees = t8_cmesh_get_num_trees (cmesh);
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
    if (num_trees > 1) {
      /* No random refinement beyond the uniform level, since the mesh has more elements */
      data.random_level = 1;
    }
    /* A point in general position inside the tree (single tree, with decreasing coordinates, which lie in the
     * reference elements of all eclasses), or near the center of the mesh, which is shared by several trees.
     * The deep elements then touch other trees across faces, edges and vertices. */
    double lower[3] = { 1e9, 1e9, 1e9 }, upper[3] = { -1e9, -1e9, -1e9 };
    for (t8_gloidx_t itree = 0; itree < (num_trees > 1 ? num_trees : 1); itree++) {
      const double *vertices = t8_cmesh_get_tree_vertices (cmesh, itree);
      for (int ivertex = 0; ivertex < t8_eclass_num_vertices[t8_cmesh_get_tree_class (cmesh, itree)]; ivertex++) {
        for (int idim = 0; idim < dim; idim++) {
          lower[idim] = SC_MIN (lower[idim], vertices[3 * ivertex + idim]);
          upper[idim] = SC_MAX (upper[idim], vertices[3 * ivertex + idim]);
        }
      }
    }
    const double single_tree_point[3] = { 0.6, 0.4, 0.2 };
    const double offsets[3] = { 0.3, 0.2, 0.1 };
    data.point = { 0, 0, 0 };
    for (int idim = 0; idim < dim; idim++) {
      const double extent = upper[idim] - lower[idim];
      if (num_trees > 1) {
        /* The offset is a fraction of the smallest element length */
        data.point[idim] = lower[idim] + extent / 2 + offsets[idim] * extent / 2 / ((int64_t) 1 << data.deep_level);
      }
      else {
        data.point[idim] = lower[idim] + single_tree_point[idim] * extent;
      }
    }
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
    const int min_dim = t8_test_min_dim (type, data.mesh.dim);
    int balanced = 1;
    for (const t8_test_leaf &leaf : t8_test_local_leaves (forest)) {
      for (const t8_test_leaf &other : all_leaves) {
        if (std::abs (leaf.level - other.level) > 1 && data.mesh.touch_dim (leaf, other) >= min_dim) {
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
  const int min_dim = t8_test_min_dim (type, data.mesh.dim);
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
      if (data.mesh.touch_dim (leaf, all_leaves[iother]) >= min_dim) {
        expected.insert (all_leaves[iother]);
        break;
      }
    }
  }

  std::set<t8_test_leaf> ghosts;
  for (t8_locidx_t itree = 0; itree < t8_forest_ghost_num_trees (unbalanced); itree++) {
    const t8_eclass_t eclass = t8_forest_ghost_get_tree_class (unbalanced, itree);
    const t8_gloidx_t gtreeid = t8_forest_ghost_get_global_treeid (unbalanced, itree);
    for (t8_locidx_t ielem = 0; ielem < t8_forest_ghost_tree_num_leaf_elements (unbalanced, itree); ielem++) {
      ghosts.insert (
        t8_test_leaf_from_element (eclass, gtreeid, t8_forest_ghost_get_leaf_element (unbalanced, itree, ielem)));
    }
  }
  EXPECT_EQ (ghosts.size (), expected.size ());
  for (const t8_test_leaf &leaf : expected) {
    EXPECT_TRUE (ghosts.count (leaf)) << "Missing ghost in tree " << leaf.gtreeid << " at level " << leaf.level
                                      << " type " << leaf.type << " anchor " << leaf.coords[0] << " "
                                      << leaf.coords[1] << " " << leaf.coords[2];
  }
  for (const t8_test_leaf &leaf : ghosts) {
    EXPECT_TRUE (expected.count (leaf)) << "Ghost that does not touch a local leaf in tree " << leaf.gtreeid
                                        << " at level " << leaf.level;
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

static t8_cmesh_t
t8_test_pad (const t8_eclass_t eclass, const int num)
{
  const double boundary[24] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 };
  const int dim = t8_eclass_to_dimension[eclass];
  return t8_cmesh_new_hypercube_pad (eclass, sc_MPI_COMM_WORLD, boundary, num, dim > 1 ? num : 0, dim > 2 ? num : 0, 0);
}

static const t8_test_cmesh t8_test_cmeshes[] = {
  { "Line", [] () { return t8_cmesh_new_from_class (T8_ECLASS_LINE, sc_MPI_COMM_WORLD); } },
  { "Quad", [] () { return t8_cmesh_new_from_class (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD); } },
  { "Triangle", [] () { return t8_cmesh_new_from_class (T8_ECLASS_TRIANGLE, sc_MPI_COMM_WORLD); } },
  { "Hex", [] () { return t8_cmesh_new_from_class (T8_ECLASS_HEX, sc_MPI_COMM_WORLD); } },
  { "Tet", [] () { return t8_cmesh_new_from_class (T8_ECLASS_TET, sc_MPI_COMM_WORLD); } },
  { "Prism", [] () { return t8_cmesh_new_from_class (T8_ECLASS_PRISM, sc_MPI_COMM_WORLD); } },
  { "Pyramid", [] () { return t8_cmesh_new_from_class (T8_ECLASS_PYRAMID, sc_MPI_COMM_WORLD); } },
  { "Lines", [] () { return t8_test_pad (T8_ECLASS_LINE, 2); } },
  { "Quads", [] () { return t8_test_pad (T8_ECLASS_QUAD, 2); } },
  { "Triangles", [] () { return t8_test_pad (T8_ECLASS_TRIANGLE, 2); } },
  { "Hexes", [] () { return t8_test_pad (T8_ECLASS_HEX, 2); } },
  { "Tets", [] () { return t8_cmesh_new_hypercube (T8_ECLASS_TET, sc_MPI_COMM_WORLD, 0, 0, 0); } },
  { "Prisms", [] () { return t8_test_pad (T8_ECLASS_PRISM, 2); } },
  { "Pyramids", [] () { return t8_cmesh_new_hypercube (T8_ECLASS_PYRAMID, sc_MPI_COMM_WORLD, 0, 0, 0); } },
  { "Hybrid", [] () { return t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0); } },
};

auto print_cmesh_ghost_type = [] (const testing::TestParamInfo<std::tuple<t8_test_cmesh, t8_ghost_type_t>> &info) {
  static const char *type_names[] = { "none", "faces", "edges", "vertices" };
  return std::get<0> (info.param).first + "_" + type_names[std::get<1> (info.param)];
};

INSTANTIATE_TEST_SUITE_P (t8_gtest_ghost_balance_boundaries, forest_ghost_balance_boundaries,
                          ::testing::Combine (::testing::ValuesIn (t8_test_cmeshes),
                                              ::testing::Values (T8_GHOST_FACES, T8_GHOST_EDGES, T8_GHOST_VERTICES)),
                          print_cmesh_ghost_type);

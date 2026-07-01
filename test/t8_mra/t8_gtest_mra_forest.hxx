#ifndef T8_GTEST_MRA_FOREST_HXX
#define T8_GTEST_MRA_FOREST_HXX

#ifdef T8_ENABLE_MRA

#include <gtest/gtest.h>

#include <t8.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_mra/t8_mra.hxx>

#include <array>
#include <cmath>
#include <string>

namespace mra_test
{

template <t8_eclass TShape, int U_, int P_>
struct Config
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int U = U_;
  static constexpr int P = P_;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;
};

using Configs
  = ::testing::Types<Config<T8_ECLASS_LINE, 1, 1>, Config<T8_ECLASS_LINE, 1, 2>, Config<T8_ECLASS_LINE, 1, 3>,
                     Config<T8_ECLASS_LINE, 1, 4>, Config<T8_ECLASS_LINE, 2, 2>, Config<T8_ECLASS_TRIANGLE, 1, 1>,
                     Config<T8_ECLASS_TRIANGLE, 1, 2>,
                     Config<T8_ECLASS_TRIANGLE, 1, 3>, Config<T8_ECLASS_TRIANGLE, 1, 4>,
                     Config<T8_ECLASS_TRIANGLE, 2, 2>, Config<T8_ECLASS_QUAD, 1, 1>, Config<T8_ECLASS_QUAD, 1, 2>,
                     Config<T8_ECLASS_QUAD, 1, 3>, Config<T8_ECLASS_QUAD, 1, 4>, Config<T8_ECLASS_QUAD, 2, 2>,
                     Config<T8_ECLASS_HEX, 1, 2>, Config<T8_ECLASS_HEX, 1, 3>, Config<T8_ECLASS_HEX, 2, 2>>;

struct ConfigNames
{
  template <typename T>
  static std::string
  GetName (int)
  {
    return std::string (t8_eclass_to_string[T::Shape]) + "_U" + std::to_string (T::U) + "_P" + std::to_string (T::P);
  }
};

/* Smooth test function */
template <int U, int DIM>
auto
smooth_func ()
{
  if constexpr (DIM == 1)
    return [] (double x) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = std::sin (2.0 * M_PI * (u + 1) * x);
      return res;
    };
  else if constexpr (DIM == 2)
    return [] (double x, double y) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = std::sin (2.0 * M_PI * (u + 1) * x) * std::sin (2.0 * M_PI * y);
      return res;
    };
  else
    return [] (double x, double y, double z) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = std::sin (2.0 * M_PI * (u + 1) * x) * std::sin (2.0 * M_PI * y) * std::sin (2.0 * M_PI * z);
      return res;
    };
}

/* Discontinuous test function: jump along a circle/sphere segment. Each side is
 * a polynomial of total degree P-1, hence exactly representable at order P */
template <int U, int P, int DIM>
auto
jump_func ()
{
  constexpr int d = P - 1;
  if constexpr (DIM == 1)
    return [] (double x) {
      std::array<double, U> res;
      const double r = x * x;
      const double in = std::pow (0.5 + x, d);
      const double out = std::pow (0.3 + x, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * ((r < 0.25) ? (3.0 + in) : out);
      return res;
    };
  else if constexpr (DIM == 2)
    return [] (double x, double y) {
      std::array<double, U> res;
      const double r = x * x + y * y;
      const double in = std::pow (0.5 + x + y, d);
      const double out = std::pow (0.3 + x - y, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * ((r < 0.25) ? (3.0 + in) : out);
      return res;
    };
  else
    return [] (double x, double y, double z) {
      std::array<double, U> res;
      const double r = x * x + y * y + z * z;
      const double in = std::pow (0.5 + x + y + z, d);
      const double out = std::pow (0.3 + x - y + z, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * ((r < 0.25) ? (3.0 + in) : out);
      return res;
    };
}

/* Constant initial data: all details vanish exactly. */
template <int U, int DIM>
auto
constant_func (double amplitude = 3.0)
{
  if constexpr (DIM == 1)
    return [amplitude] (double) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = amplitude * (u + 1);
      return res;
    };
  else if constexpr (DIM == 2)
    return [amplitude] (double, double) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = amplitude * (u + 1);
      return res;
    };
  else
    return [amplitude] (double, double, double) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = amplitude * (u + 1);
      return res;
    };
}

/* Single-piece polynomial of total degree P-1, hence exactly representable in
 * the order-P space on every shape (Dubiner: total degree <= P-1; tensor
 * Legendre: per-variable degree <= P-1). No jump, so a correct projection
 * leaves zero details. */
template <int U, int P, int DIM>
auto
poly_func ()
{
  constexpr int d = P - 1;
  if constexpr (DIM == 1)
    return [] (double x) {
      std::array<double, U> res;
      const double base = std::pow (0.5 + 0.3 * x, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * base;
      return res;
    };
  else if constexpr (DIM == 2)
    return [] (double x, double y) {
      std::array<double, U> res;
      const double base = std::pow (0.5 + 0.3 * x + 0.4 * y, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * base;
      return res;
    };
  else
    return [] (double x, double y, double z) {
      std::array<double, U> res;
      const double base = std::pow (0.5 + 0.3 * x + 0.4 * y + 0.2 * z, d);
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * base;
      return res;
    };
}

/* Every forest leaf must have its lmi stored with the element's level, and
 * exactly the leaves must be present in lmi_map. */
template <typename MRA>
void
expect_forest_map_consistent (MRA &mra)
{
  auto *forest = mra.get_forest ();
  auto *user_data = mra.get_user_data ();
  auto *lmi_map = mra.get_lmi_map ();
  const auto *scheme = t8_forest_get_scheme (forest);

  ASSERT_EQ (static_cast<size_t> (t8_forest_get_local_num_leaf_elements (forest)), lmi_map->size ());

  auto current_idx = t8_locidx_t { 0 };
  const auto num_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
    const auto tree_class = t8_forest_get_tree_class (forest, tree_idx);
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree_idx);

    for (t8_locidx_t ele_idx = 0; ele_idx < num_elements; ++ele_idx, ++current_idx) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree_idx, ele_idx);
      const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, current_idx);

      EXPECT_EQ (lmi.level (), static_cast<unsigned int> (scheme->element_get_level (tree_class, element)));
      EXPECT_TRUE (lmi_map->contains (lmi));
    }
  }
}

template <typename MapT>
void
expect_maps_equal (const MapT &expected, const MapT &actual, unsigned int max_level, double tol)
{
  for (auto l = 0u; l <= max_level; ++l) {
    ASSERT_EQ (expected[l].size (), actual[l].size ()) << "entry count differs on level " << l;

    for (const auto &[lmi, data] : expected[l]) {
      ASSERT_TRUE (actual.contains (lmi)) << "missing lmi on level " << l;

      const auto &other = actual.get (lmi);
      ASSERT_EQ (data.u_coeffs.size (), other.u_coeffs.size ());
      for (auto i = 0u; i < data.u_coeffs.size (); ++i)
        EXPECT_NEAR (data.u_coeffs[i], other.u_coeffs[i], tol) << "coeff " << i << " on level " << l;
    }
  }
}

/* (cmesh, scheme, multiscale) triple on a unit hypercube */
template <t8_eclass TShape, int U, int P>
class mra_example {
 public:
  using multiscale = t8_mra::multiscale<TShape, U, P>;
  using levelmultiindex = typename multiscale::levelmultiindex;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;

  explicit mra_example (int max_level): mra (max_level, sc_MPI_COMM_WORLD), max_level (max_level)
  {
    cmesh = t8_cmesh_new_hypercube (TShape, sc_MPI_COMM_WORLD, 0, 0, 0);
    scheme = t8_scheme_new_default ();
    t8_cmesh_ref (cmesh);
    t8_scheme_ref (const_cast<t8_scheme *> (scheme));
  }

  ~mra_example ()
  {
    mra.cleanup ();
    t8_cmesh_destroy (&cmesh);
    t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
  }

  mra_example (const mra_example &) = delete;

  mra_example &
  operator= (const mra_example &) = delete;

  template <typename F>
  void
  init (F &&f)
  {
    mra.initialize_data (cmesh, scheme, max_level, std::forward<F> (f));
  }

  template <typename F>
  void
  init_adaptive (F &&f)
  {
    mra.initialize_data_adaptive (cmesh, scheme, max_level, std::forward<F> (f));
  }

  multiscale *
  operator->()
  {
    return &mra;
  }

  multiscale &
  operator* ()
  {
    return mra;
  }

  multiscale mra;
  t8_cmesh_t cmesh;
  const t8_scheme *scheme;
  int max_level;
};

/* Grading invariant: across every face the leaf levels differ by at most
 * `slack` (slack 1 = the 2:1 graded grid). Domain-boundary faces carry no
 * neighbour and are skipped. Serial only -- a process boundary without a ghost
 * layer would also report zero neighbours. */
template <typename MRA>
void
expect_grid_graded (MRA &mra, int slack = 1)
{
  auto *forest = mra.get_forest ();
  const auto *scheme = t8_forest_get_scheme (forest);
  const auto num_trees = t8_forest_get_num_local_trees (forest);

  for (t8_locidx_t tree = 0; tree < num_trees; ++tree) {
    const auto tree_class = t8_forest_get_tree_class (forest, tree);
    const auto num_elements = t8_forest_get_tree_num_leaf_elements (forest, tree);

    for (t8_locidx_t e = 0; e < num_elements; ++e) {
      const auto *element = t8_forest_get_leaf_element_in_tree (forest, tree, e);
      const int level = scheme->element_get_level (tree_class, element);
      const int num_faces = scheme->element_get_num_faces (tree_class, element);

      for (int face = 0; face < num_faces; ++face) {
        const t8_element_t **neighbors = nullptr;
        int *dual_faces = nullptr;
        t8_locidx_t *element_indices = nullptr;
        int num_neighbors = 0;
        t8_eclass_t neigh_class;

        t8_forest_leaf_face_neighbors (forest, tree, element, &neighbors, face, &dual_faces, &num_neighbors,
                                       &element_indices, &neigh_class);

        for (int n = 0; n < num_neighbors; ++n) {
          const int neigh_level = scheme->element_get_level (neigh_class, neighbors[n]);
          EXPECT_LE (std::abs (level - neigh_level), slack) << "ungraded face neighbour at level " << level;
        }

        if (num_neighbors > 0) {
          T8_FREE (neighbors);
          T8_FREE (element_indices);
          T8_FREE (dual_faces);
        }
      }
    }
  }
}

}  // namespace mra_test

#endif /* T8_ENABLE_MRA */

#endif /* T8_GTEST_MRA_FOREST_HXX */

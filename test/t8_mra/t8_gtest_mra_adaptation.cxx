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

/** \file t8_gtest_mra_adaptation.cxx
 * Regression tests for the MRA path (multiscale_base + mst +
 * multiscale_adaptation), typed over (element shape, U components, order P):
 *   - forward/inverse multiscale transformation round-trip is the identity
 *   - coarsen -> refine -> coarsen returns to the coarsened grid and data
 *   - lmi_map mirrors the forest leaves exactly after every adaptation step
 *
 * Shapes: TRIANGLE (hardcoded masks, P = 1..4), QUAD and HEX (computed
 * masks, arbitrary P).
 */

#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#include <t8_mra/t8_mra.hxx>

#include <array>
#include <cmath>
#include <string>

namespace
{

template <t8_eclass TShape, int U_, int P_>
struct Config
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int U = U_;
  static constexpr int P = P_;
  static constexpr int DIM = t8_mra::t8_eclass_dim<TShape> ();
};

/* Smooth test function for the MST round-trip; components differ to catch
 * component-indexing bugs. */
template <int U, int DIM>
auto
smooth_func ()
{
  if constexpr (DIM == 2)
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

/* Discontinuous test function: jump along a circle/sphere segment.
 * Guarantees significant details at the jump, hence coarsening keeps a fine
 * band there and refinement must trigger via the neighbour prediction. */
template <int U, int DIM>
auto
jump_func ()
{
  if constexpr (DIM == 2)
    return [] (double x, double y) {
      std::array<double, U> res;
      const double r = x * x + y * y;
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * ((r < 0.25) ? (x * y + x + 3.0) : (x * x * y - 2.0 * x * y * y + 3.0 * x));
      return res;
    };
  else
    return [] (double x, double y, double z) {
      std::array<double, U> res;
      const double r = x * x + y * y + z * z;
      for (auto u = 0; u < U; ++u)
        res[u] = (u + 1) * ((r < 0.25) ? (x * y + z + 3.0) : (x * x * y - 2.0 * y * z + 3.0 * x));
      return res;
    };
}

template <t8_eclass TShape, int U, int P>
t8_mra::multiscale<TShape, U, P>
make_mra (int max_level)
{
  return t8_mra::multiscale<TShape, U, P> (max_level, sc_MPI_COMM_WORLD);
}

/* Every forest leaf must have its lmi stored in lmi_idx with the element's
 * level, and exactly the leaves must be present in lmi_map. */
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

/* Constant initial data: all detail coefficients vanish exactly, so the
 * bottom-up initialization must never keep a refined level. */
template <int U, int DIM>
auto
constant_func ()
{
  if constexpr (DIM == 2)
    return [] (double, double) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = 3.0 * (u + 1);
      return res;
    };
  else
    return [] (double, double, double) {
      std::array<double, U> res;
      for (auto u = 0; u < U; ++u)
        res[u] = 3.0 * (u + 1);
      return res;
    };
}

/* Custom coarsening criterion without prepare(): nothing is significant
 * -> coarsening must collapse the grid completely. Exercises the
 * coarsening_criterion extension point. */
struct collapse_criterion
{
  template <typename MRA>
  bool
  significant (MRA &, const typename MRA::levelmultiindex &)
  {
    return false;
  }
};

template <typename Cfg>
class mra_adaptation: public ::testing::Test {
};

using Configs = ::testing::Types<
  /* Triangle: hardcoded masks, P = 1..4 */
  Config<T8_ECLASS_TRIANGLE, 1, 1>, Config<T8_ECLASS_TRIANGLE, 1, 2>, Config<T8_ECLASS_TRIANGLE, 1, 3>,
  Config<T8_ECLASS_TRIANGLE, 1, 4>, Config<T8_ECLASS_TRIANGLE, 2, 2>,
  /* Quad: computed masks */
  Config<T8_ECLASS_QUAD, 1, 1>, Config<T8_ECLASS_QUAD, 1, 2>, Config<T8_ECLASS_QUAD, 1, 3>,
  Config<T8_ECLASS_QUAD, 1, 4>, Config<T8_ECLASS_QUAD, 2, 2>,
  /* Hex (3D): computed masks */
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

TYPED_TEST_SUITE (mra_adaptation, Configs, ConfigNames);

/* Forward then inverse MST over the full level range must reproduce the
 * single-scale data up to round-off. */
TYPED_TEST (mra_adaptation, mst_roundtrip)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;
  auto mra = make_mra<Shape, U, P> (max_level);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  /* The forest takes ownership of cmesh and scheme; keep our own references */
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data (cmesh, scheme, max_level, smooth_func<U, DIM> ());
  expect_forest_map_consistent (mra);

  const auto snapshot = *mra.get_lmi_map ();

  mra.multiscale_transformation (0, max_level);
  /* All single-scale data restricted to level 0 */
  EXPECT_EQ (mra.get_lmi_map ()->size (), mra.get_lmi_map ()->size (0));

  mra.inverse_multiscale_transformation (0, max_level);
  expect_maps_equal (snapshot, *mra.get_lmi_map (), max_level, 1e-9);

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
}

/* Coarsening discards only non-significant details; refinement adds only
 * zero-detail children (plus Harten grading). A second coarsening therefore
 * must return exactly to the first coarsened grid and data. */
TYPED_TEST (mra_adaptation, coarsen_refine_roundtrip)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  if constexpr (P == 1) {
    /* With piecewise constants the approximation error of the smooth pieces
     * dominates the details on all levels; whether coarsening triggers is a
     * threshold tuning question, not an algorithmic invariant. */
    GTEST_SKIP () << "coarsening behaviour for P=1 is threshold-dependent";
  }
  else {
    const int max_level = (DIM == 3) ? 4 : 5;
    auto mra = make_mra<Shape, U, P> (max_level);

    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
    auto *scheme = t8_scheme_new_default ();
    t8_cmesh_ref (cmesh);
    t8_scheme_ref (const_cast<t8_scheme *> (scheme));

    mra.initialize_data (cmesh, scheme, max_level, jump_func<U, DIM> ());
    const auto num_uniform = t8_forest_get_global_num_leaf_elements (mra.get_forest ());

    mra.coarsen (0, max_level);
    expect_forest_map_consistent (mra);
    const auto num_coarse = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    EXPECT_LT (num_coarse, num_uniform) << "smooth regions must coarsen";

    const auto snapshot = *mra.get_lmi_map ();

    mra.refine (0, max_level);
    expect_forest_map_consistent (mra);
    const auto num_refined = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    EXPECT_GT (num_refined, num_coarse) << "neighbour prediction must refine around the jump";

    mra.coarsen (0, max_level);
    expect_forest_map_consistent (mra);
    const auto num_recoarse = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    EXPECT_EQ (num_recoarse, num_coarse) << "zero-detail children must coarsen away again";

    expect_maps_equal (snapshot, *mra.get_lmi_map (), max_level, 1e-8);

    mra.cleanup ();
    t8_cmesh_destroy (&cmesh);
    t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
  }
}

/* A criterion that finds nothing significant must collapse the grid to the
 * base cells, regardless of the data. */
TYPED_TEST (mra_adaptation, custom_criterion_collapse)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;
  auto mra = make_mra<Shape, U, P> (max_level);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data (cmesh, scheme, max_level, jump_func<U, DIM> ());

  mra.coarsen (0, max_level, collapse_criterion {});
  expect_forest_map_consistent (mra);

  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (mra.get_forest ()),
             t8_forest_get_num_global_trees (mra.get_forest ()));

  mra.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
}

/* Bottom-up initialization on constant data: every detail is exactly zero,
 * so nothing is significant and the result must be the level-1 grid with
 * the level-1 projection. */
TYPED_TEST (mra_adaptation, bottom_up_constant_collapses)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;
  auto mra = make_mra<Shape, U, P> (max_level);
  auto reference = make_mra<Shape, U, P> (max_level);

  t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
  auto *scheme = t8_scheme_new_default ();
  /* Two forests consume one reference each; keep our own on top */
  t8_cmesh_ref (cmesh);
  t8_cmesh_ref (cmesh);
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));
  t8_scheme_ref (const_cast<t8_scheme *> (scheme));

  mra.initialize_data_adaptive (cmesh, scheme, max_level, constant_func<U, DIM> ());
  expect_forest_map_consistent (mra);

  reference.initialize_data (cmesh, scheme, 1, constant_func<U, DIM> ());

  EXPECT_EQ (t8_forest_get_global_num_leaf_elements (mra.get_forest ()),
             t8_forest_get_global_num_leaf_elements (reference.get_forest ()));
  expect_maps_equal (*reference.get_lmi_map (), *mra.get_lmi_map (), max_level, 1e-9);

  mra.cleanup ();
  reference.cleanup ();
  t8_cmesh_destroy (&cmesh);
  t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
}

/* Bottom-up initialization on discontinuous data must produce an adaptive
 * grid: strictly coarser than uniform, refined to max_level at the jump,
 * and a valid starting point for the regular adaptation cycle. */
TYPED_TEST (mra_adaptation, bottom_up_adaptive_grid)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  if constexpr (P == 1) {
    GTEST_SKIP () << "coarsening behaviour for P=1 is threshold-dependent";
  }
  else {
    using levelmultiindex = typename t8_mra::multiscale<Shape, U, P>::levelmultiindex;

    const int max_level = (DIM == 3) ? 4 : 5;
    auto mra = make_mra<Shape, U, P> (max_level);

    t8_cmesh_t cmesh = t8_cmesh_new_hypercube (Shape, sc_MPI_COMM_WORLD, 0, 0, 0);
    auto *scheme = t8_scheme_new_default ();
    t8_cmesh_ref (cmesh);
    t8_scheme_ref (const_cast<t8_scheme *> (scheme));

    mra.initialize_data_adaptive (cmesh, scheme, max_level, jump_func<U, DIM> ());
    expect_forest_map_consistent (mra);

    const auto num_adaptive = t8_forest_get_global_num_leaf_elements (mra.get_forest ());
    const auto num_trees = t8_forest_get_num_global_trees (mra.get_forest ());
    const auto num_uniform
      = num_trees * static_cast<t8_gloidx_t> (std::pow (levelmultiindex::NUM_CHILDREN, max_level));

    EXPECT_LT (num_adaptive, num_uniform) << "smooth regions must stay coarse";
    EXPECT_GT (mra.get_lmi_map ()->size (max_level), 0u) << "the jump must reach max_level";
    EXPECT_EQ (mra.get_lmi_map ()->size (0), 0u) << "level 1 is the minimum level";

    /* The result must be a valid input for the regular adaptation cycle */
    mra.coarsen (1, max_level);
    expect_forest_map_consistent (mra);
    EXPECT_LE (t8_forest_get_global_num_leaf_elements (mra.get_forest ()), num_adaptive);

    mra.cleanup ();
    t8_cmesh_destroy (&cmesh);
    t8_scheme_unref (const_cast<t8_scheme **> (&scheme));
  }
}

}  // namespace

#endif /* T8_ENABLE_MRA */

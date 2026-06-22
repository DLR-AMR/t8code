#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include "t8_gtest_mra_forest.hxx"

namespace
{

using namespace mra_test;

constexpr double eps = 1e-12;

template <typename Cfg>
class mra_criteria: public ::testing::Test {};

TYPED_TEST_SUITE (mra_criteria, Configs, ConfigNames);

/* Uniform max-level grid with computed details: the common starting point for
 * every threshold check. */
template <typename Case, typename F>
void
init_and_decompose (Case &c, F &&f)
{
  c.init (std::forward<F> (f));
  c->multiscale_decomposition (0, c.max_level);
}

/* threshold_scaling_factor is a domain integral clamped to >= 1, and the DG
 * projection is linear, so doubling the data doubles the (unclamped) factor. */
TYPED_TEST (mra_criteria, threshold_scaling_factor_is_clamped_and_linear)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 2 : 3;

  /* Vanishing data -> the clamp pins every component at 1. */
  {
    mra_example<Shape, U, P> example (max_level);
    example.init (constant_func<U, DIM> (1e-9));
    const auto factor = example->threshold_scaling_factor ();
    for (auto u = 0u; u < U; ++u)
      EXPECT_NEAR (factor[u], 1.0, eps) << "tiny data must clamp to 1, component " << u;
  }

  /* Large data -> unclamped, and linear in the amplitude. */
  mra_example<Shape, U, P> c1 (max_level);
  mra_example<Shape, U, P> c2 (max_level);
  c1.init (constant_func<U, DIM> (100.0));
  c2.init (constant_func<U, DIM> (200.0));

  const auto f1 = c1->threshold_scaling_factor ();
  const auto f2 = c2->threshold_scaling_factor ();
  for (auto u = 0u; u < U; ++u) {
    ASSERT_GT (f1[u], 1.0) << "amplitude 100 should exceed the clamp, component " << u;
    EXPECT_NEAR (f2[u] / f1[u], 2.0, eps) << "factor must scale linearly, component " << u;
  }
}

/* local_threshold_value depends only on (vol, level). */
TYPED_TEST (mra_criteria, local_threshold_value_scales_as_sqrt_num_children)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using LMI = typename t8_mra::multiscale<Shape, U, P>::levelmultiindex;

  const int max_level = (DIM == 3) ? 3 : 4;
  const double ratio = std::sqrt (static_cast<double> (LMI::NUM_CHILDREN));

  mra_example<Shape, U, P> example (max_level);
  init_and_decompose (example, jump_func<U, P, DIM> ());
  auto &mra = example.mra;

  std::size_t pairs = 0;
  for (const int gamma : { 1, 2, 3 })
    for (auto l = 2u; l <= static_cast<unsigned int> (max_level); ++l)
      for (const auto &[lmi, detail] : mra.d_map[l]) {
        const auto par = LMI::parent (lmi);
        if (!mra.d_map.contains (par))
          continue;
        ++pairs;
        EXPECT_NEAR (mra.local_threshold_value (lmi, gamma), ratio * mra.local_threshold_value (par, gamma),
                     eps * mra.local_threshold_value (par, gamma))
          << "level " << l << " gamma " << gamma;
      }
  EXPECT_GT (pairs, 0u) << "decomposition must leave parent/child detail pairs to compare";
}

/* scaled_detail_norm divides the raw detail norm by c_scaling componentwise */
TYPED_TEST (mra_criteria, scaled_detail_norm_respects_c_scaling)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> example (max_level);
  init_and_decompose (example, jump_func<U, P, DIM> ());
  auto &mra = example.mra;

  bool has_nonzero = false;
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, detail] : mra.d_map[l]) {
      mra.c_scaling.fill (1.0);
      const double base = mra.scaled_detail_norm (lmi);
      if (base <= eps)
        continue;
      has_nonzero = true;

      mra.c_scaling.fill (4.0);
      EXPECT_NEAR (mra.scaled_detail_norm (lmi), base / 4.0, eps * base) << "level " << l;
    }
  EXPECT_TRUE (has_nonzero) << "the jump must leave at least one nonzero detail";
}

/* Raising the threshold constant can only remove leaves from the significant
 * set. */
TYPED_TEST (mra_criteria, hard_thresholding_significant_set_shrinks_with_threshold)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> example (max_level);
  example.init (jump_func<U, P, DIM> ());
  auto &mra = example.mra;

  t8_mra::hard_thresholding low { 0.1, 2 };
  t8_mra::hard_thresholding high { 10.0, 2 };
  low.prepare (mra);
  high.prepare (mra);

  mra.multiscale_decomposition (0, max_level);

  std::size_t number_low = 0, number_high = 0;
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, detail] : mra.d_map[l]) {
      const bool significant_low = low.significant (mra, lmi);
      const bool significant_high = high.significant (mra, lmi);

      number_low += significant_low;
      number_high += significant_high;
      EXPECT_TRUE (!significant_high || significant_low)
        << "a higher threshold must not gain significance, level " << l;
    }
  EXPECT_LE (number_high, number_low);
  EXPECT_GT (number_low, 0u) << "the jump must make some family significant at a low threshold";
}

/* The Harten steep-gradient refinement test is 2^(P+1) times stricter than the
 * neighbour-grading test. */
TYPED_TEST (mra_criteria, harten_refine_implies_refine_neighbours)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> example (max_level);
  example.init (jump_func<U, P, DIM> ());
  auto &mra = example.mra;

  t8_mra::harten_prediction crit { 1.0, 2 };
  crit.prepare (mra);

  mra.multiscale_decomposition (0, max_level);

  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, detail] : mra.d_map[l])
      if (crit.refine (mra, lmi))
        EXPECT_TRUE (crit.refine_neighbours (mra, lmi)) << "refine must imply grading, level " << l;
}

}  // namespace

#endif /* T8_ENABLE_MRA */

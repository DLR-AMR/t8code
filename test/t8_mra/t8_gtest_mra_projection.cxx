#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include "t8_gtest_mra_forest.hxx"

namespace
{

using namespace mra_test;

constexpr double eps = 1e-12;

template <typename Cfg>
class mra_projection: public ::testing::Test {};

TYPED_TEST_SUITE (mra_projection, Configs, ConfigNames);

/* Constant data has no variation, so every non-constant DG mode must vanish on
 * every leaf. */
TYPED_TEST (mra_projection, constant_projects_to_pure_constant_mode)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using element_t = typename t8_mra::multiscale<Shape, U, P>::element_t;

  const int max_level = (DIM == 3) ? 2 : 3;

  mra_example<Shape, U, P> example (max_level);
  example.init (constant_func<U, DIM> (2.5));
  auto &mra = example.mra;

  auto *lmi_map = mra.get_lmi_map ();
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, data] : (*lmi_map)[l])
      for (auto u = 0u; u < U; ++u)
        for (auto i = 1u; i < element_t::DOF; ++i)
          EXPECT_NEAR (data.u_coeffs[element_t::dg_idx (u, i)], 0.0, eps)
            << "non-constant mode " << i << " must vanish for constant data, component " << u;
}

/* A polynomial of total degree <= P-1 lives in the coarse space, so projecting
 * it on a uniform grid and running a full multiscale decomposition leaves zero
 * details -- the end-to-end projection exactness statement (mirrors
 * multilaepsch's project / cancellation, but through the production forest
 * projection path and, on triangles, the derived vertex ordering). */
TYPED_TEST (mra_projection, representable_polynomial_has_no_details)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using detail_t = typename t8_mra::multiscale<Shape, U, P>::detail_t;

  const int max_level = (DIM == 3) ? 3 : 4;

  mra_example<Shape, U, P> example (max_level);
  example.init (poly_func<U, P, DIM> ());
  auto &mra = example.mra;

  mra.multiscale_decomposition (0, max_level);

  std::size_t families = 0;
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, detail] : mra.d_map[l]) {
      ++families;
      for (const double d : detail.d_coeffs)
        EXPECT_NEAR (d, 0.0, eps) << "representable polynomial must leave no detail, level " << l;
    }
  EXPECT_GT (families, 0u) << "decomposition must produce families to check";
}

/* The DG projection is a linear operator, so projecting the scaled data yields
 * the scaled coefficients on every leaf. */
TYPED_TEST (mra_projection, projection_is_linear)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using element_t = typename t8_mra::multiscale<Shape, U, P>::element_t;

  const int max_level = (DIM == 3) ? 2 : 3;
  constexpr double scale = 3.0;

  auto base = jump_func<U, P, DIM> ();

  mra_example<Shape, U, P> plain (max_level);
  mra_example<Shape, U, P> scaled (max_level);
  plain.init (base);
  scaled.init ([base] (auto... x) {
    auto v = base (x...);
    for (auto &e : v)
      e *= scale;
    return v;
  });

  auto *plain_map = plain->get_lmi_map ();
  auto *scaled_map = scaled->get_lmi_map ();
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, data] : (*plain_map)[l]) {
      ASSERT_TRUE (scaled_map->contains (lmi)) << "both grids must be identical, level " << l;
      const auto &other = scaled_map->get (lmi);
      for (auto u = 0u; u < U; ++u)
        for (auto i = 0u; i < element_t::DOF; ++i) {
          const auto idx = element_t::dg_idx (u, i);
          EXPECT_NEAR (other.u_coeffs[idx], scale * data.u_coeffs[idx], eps + eps * std::abs (data.u_coeffs[idx]))
            << "coeff " << i << " component " << u << " level " << l;
        }
    }
}

}  // namespace

#endif /* T8_ENABLE_MRA */

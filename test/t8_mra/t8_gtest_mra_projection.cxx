#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include "t8_gtest_mra_forest.hxx"

#include <t8_forest/t8_forest_geometrical.h>

#include <set>
#include <vector>

namespace
{

using namespace mra_test;

constexpr double eps = 1e-12;

/* Sample points guaranteed inside (or on) a leaf: every vertex plus the
 * centroid. A degree <= P-1 polynomial is reproduced exactly at all of them. */
template <int DIM>
std::vector<std::array<double, DIM>>
leaf_sample_points (t8_forest_t forest, t8_locidx_t tree_idx, const t8_element_t *element, int num_vertices)
{
  std::vector<std::array<double, DIM>> samples;
  std::array<double, DIM> centroid = {};

  for (int v = 0; v < num_vertices; ++v) {
    double coord[3] = {};
    t8_forest_element_coordinate (forest, tree_idx, element, v, coord);

    std::array<double, DIM> point;
    for (auto d = 0; d < DIM; ++d) {
      point[d] = coord[d];
      centroid[d] += coord[d] / num_vertices;
    }
    samples.push_back (point);
  }
  samples.push_back (centroid);

  return samples;
}

/* Call an (x,y[,z]) test function with a point stored as an array. */
template <int DIM, typename F>
auto
eval_func (F &&f, const std::array<double, DIM> &x)
{
  if constexpr (DIM == 2)
    return f (x[0], x[1]);
  else
    return f (x[0], x[1], x[2]);
}

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

/* The constant DG mode carries the cell integral: exactly vol*u0 on cartesian
 * cells and sqrt(vol)*u0 on triangles (the orthonormal Dubiner basis scales by
 * sqrt(1/(2*vol)), so u0 = mean*sqrt(vol)). */
template <t8_eclass Shape>
double
cell_integral_factor (double vol)
{
  return (Shape == T8_ECLASS_TRIANGLE) ? std::sqrt (vol) : vol;
}

/* Summing the constant modes recovers the exact domain integral: for constant
 * data c the reconstructed mass is c times the total volume. Also pins the
 * shape-dependent constant-mode normalization (a wrong triangle factor fails
 * here). */
TYPED_TEST (mra_projection, constant_mode_reconstructs_domain_mass)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using element_t = typename t8_mra::multiscale<Shape, U, P>::element_t;

  const int max_level = (DIM == 3) ? 2 : 3;
  constexpr double amplitude = 2.5;

  mra_example<Shape, U, P> example (max_level);
  example.init (constant_func<U, DIM> (amplitude));
  auto &mra = example.mra;

  std::array<double, U> mass = {};
  double total_vol = 0.0;
  auto *lmi_map = mra.get_lmi_map ();
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, data] : (*lmi_map)[l]) {
      total_vol += data.vol;
      const auto factor = cell_integral_factor<Shape> (data.vol);
      for (auto u = 0u; u < U; ++u)
        mass[u] += factor * data.u_coeffs[element_t::dg_idx (u, 0)];
    }

  for (auto u = 0u; u < U; ++u) {
    const double expected = amplitude * (u + 1) * total_vol;
    EXPECT_NEAR (mass[u], expected, eps) << "component " << u;
  }
}

/* Projection conserves the domain integral independent of resolution: a
 * representable field projected on level L and on level L+1 reconstructs the
 * same total mass (conservation of a varying field, no analytic integral
 * needed). */
TYPED_TEST (mra_projection, projection_conserves_mass_across_levels)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;
  using element_t = typename t8_mra::multiscale<Shape, U, P>::element_t;

  const int coarse_level = (DIM == 3) ? 2 : 3;

  auto reconstructed_mass = [] (auto &mra, int max_level) {
    std::array<double, U> mass = {};
    auto *lmi_map = mra.get_lmi_map ();
    for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
      for (const auto &[lmi, data] : (*lmi_map)[l]) {
        const auto factor = cell_integral_factor<Shape> (data.vol);
        for (auto u = 0u; u < U; ++u)
          mass[u] += factor * data.u_coeffs[element_t::dg_idx (u, 0)];
      }
    return mass;
  };

  auto f = poly_func<U, P, DIM> ();
  mra_example<Shape, U, P> coarse (coarse_level);
  mra_example<Shape, U, P> fine (coarse_level + 1);
  coarse.init (f);
  fine.init (f);

  const auto mass_coarse = reconstructed_mass (coarse.mra, coarse_level);
  const auto mass_fine = reconstructed_mass (fine.mra, coarse_level + 1);
  for (auto u = 0u; u < U; ++u)
    EXPECT_NEAR (mass_fine[u], mass_coarse[u], eps) << "component " << u;
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
          EXPECT_NEAR (other.u_coeffs[idx], scale * data.u_coeffs[idx], eps)
            << "coeff " << i << " component " << u << " level " << l;
        }
    }
}

/* End-to-end projection exactness, pointwise: a polynomial of total degree
 * <= P-1 is reconstructed exactly by evaluate() at every leaf vertex and
 * centroid (mirrors multilaepsch's project/val, now that a solution-evaluation
 * API exists). */
TYPED_TEST (mra_projection, projection_reconstructs_field_pointwise)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 2 : 3;
  const int num_vertices = t8_eclass_num_vertices[Shape];

  auto f = poly_func<U, P, DIM> ();
  mra_example<Shape, U, P> example (max_level);
  example.init (f);
  auto &mra = example.mra;

  auto *forest = mra.get_forest ();
  auto *user_data = mra.get_user_data ();
  auto *lmi_map = mra.get_lmi_map ();

  std::size_t checked = 0;
  mra.for_each_local_leaf (
    [&] (t8_locidx_t tree_idx, const t8_element_t *element, unsigned int local_idx, t8_gloidx_t) {
      const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_idx);
      const auto &data = lmi_map->get (lmi);

      for (const auto &point : leaf_sample_points<DIM> (forest, tree_idx, element, num_vertices)) {
        const auto got = mra.evaluate (tree_idx, element, data, point);
        const auto exact = eval_func<DIM> (f, point);
        for (auto u = 0u; u < U; ++u)
          EXPECT_NEAR (got[u], exact[u], eps) << "reconstruction must match, component " << u;
      }
      ++checked;
    });
  EXPECT_GT (checked, 0u) << "must reconstruct at least one leaf";
}

/* The triangle projection is independent of the derived vertex ordering: a
 * refined grid carries several distinct orders (ancestor.type driven), yet
 * evaluate() reproduces a representable polynomial exactly on every leaf. */
TYPED_TEST (mra_projection, triangle_projection_is_vertex_order_invariant)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  if constexpr (Shape != T8_ECLASS_TRIANGLE) {
    GTEST_SKIP () << "vertex ordering is triangle-specific";
  }
  else {
    const int max_level = 3;

    auto f = poly_func<U, P, DIM> ();
    mra_example<Shape, U, P> example (max_level);
    example.init (f);
    auto &mra = example.mra;

    auto *forest = mra.get_forest ();
    auto *user_data = mra.get_user_data ();
    auto *lmi_map = mra.get_lmi_map ();
    const int num_vertices = t8_eclass_num_vertices[Shape];

    std::set<std::array<int, 3>> orders;
    mra.for_each_local_leaf (
      [&] (t8_locidx_t tree_idx, const t8_element_t *element, unsigned int local_idx, t8_gloidx_t) {
        const auto lmi = t8_mra::get_lmi_from_forest_data (user_data, local_idx);
        const auto &data = lmi_map->get (lmi);
        orders.insert (data.order);

        for (const auto &point : leaf_sample_points<DIM> (forest, tree_idx, element, num_vertices)) {
          const auto got = mra.evaluate (tree_idx, element, data, point);
          const auto exact = eval_func<DIM> (f, point);
          for (auto u = 0u; u < U; ++u)
            EXPECT_NEAR (got[u], exact[u], eps) << "reconstruction must match under every order, component " << u;
        }
      });
    EXPECT_GT (orders.size (), 1u) << "the refined grid must exercise more than one vertex order";
  }
}

}  // namespace

#endif /* T8_ENABLE_MRA */

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

/* Constant data: only the constant mode survives on every leaf. */
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

/* Degree <= P-1 polynomial lives in the coarse space: a full decomposition
 * leaves zero details (projection exactness through the forest path; on
 * triangles also the derived vertex ordering). */
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

/* Cell integral from the constant mode: vol*u0 (cartesian), sqrt(vol)*u0
 * (triangle: Dubiner scales by sqrt(1/(2*vol))). */
template <t8_eclass Shape>
double
cell_integral_factor (double vol)
{
  return (Shape == T8_ECLASS_TRIANGLE) ? std::sqrt (vol) : vol;
}

/* Constant modes sum to the exact domain integral (c * total volume); pins the
 * shape-dependent constant-mode factor. */
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

/* Reconstructed mass is resolution-independent: same total on level L and L+1. */
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

/* Projection is linear: scaling the data scales every coefficient. */
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

/* Degree <= P-1 polynomial reconstructs exactly at every leaf vertex and
 * centroid. */
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

/* Triangle projection is independent of the derived vertex order: a refined grid
 * carries several orders, each reconstructs the polynomial exactly. */
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

/* evaluate_point finds the owning leaf: exact at interior points, nullopt
 * outside the domain. */
TYPED_TEST (mra_projection, evaluate_point_reconstructs_field)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 2 : 3;

  auto f = poly_func<U, P, DIM> ();
  mra_example<Shape, U, P> example (max_level);
  example.init (f);
  auto &mra = example.mra;

  std::vector<std::array<double, DIM>> points;
  if constexpr (DIM == 2)
    points = { { 0.3, 0.2 }, { 0.15, 0.8 }, { 0.7, 0.25 } };
  else
    points = { { 0.3, 0.2, 0.4 }, { 0.15, 0.6, 0.7 }, { 0.6, 0.25, 0.1 } };

  for (const auto &point : points) {
    const auto got = mra.evaluate_point (point);
    ASSERT_TRUE (got.has_value ()) << "a point inside the domain must be owned by a leaf";
    const auto exact = eval_func<DIM> (f, point);
    for (auto u = 0u; u < U; ++u)
      EXPECT_NEAR ((*got)[u], exact[u], eps) << "reconstruction must match, component " << u;
  }

  std::array<double, DIM> outside;
  outside.fill (0.5);
  outside[0] = 1.5;
  EXPECT_FALSE (mra.evaluate_point (outside).has_value ()) << "a point outside the domain must not be owned";
}

/* mean_val returns the cell average: for constant data it is the amplitude on
 * every leaf. */
TYPED_TEST (mra_projection, mean_val_equals_constant_data)
{
  constexpr auto Shape = TypeParam::Shape;
  constexpr auto U = TypeParam::U;
  constexpr auto P = TypeParam::P;
  constexpr auto DIM = TypeParam::DIM;

  const int max_level = (DIM == 3) ? 2 : 3;
  constexpr double amplitude = 2.5;

  mra_example<Shape, U, P> example (max_level);
  example.init (constant_func<U, DIM> (amplitude));
  auto &mra = example.mra;

  auto *lmi_map = mra.get_lmi_map ();
  for (auto l = 0u; l <= static_cast<unsigned int> (max_level); ++l)
    for (const auto &[lmi, data] : (*lmi_map)[l]) {
      const auto mean = mra.mean_val (data);
      for (auto u = 0u; u < U; ++u)
        EXPECT_NEAR (mean[u], amplitude * (u + 1), eps) << "component " << u;
    }
}

}  // namespace

#endif /* T8_ENABLE_MRA */

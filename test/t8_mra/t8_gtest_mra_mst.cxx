#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/core/mst.hxx>
#include <t8_mra/core/shape_traits.hxx>
#include <t8_mra/data/element_data.hxx>
#include <t8_mra/data/levelindex_map.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/mask_coefficients.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace
{

constexpr double eps = 1e-12;

template <t8_eclass TShape, int P_>
struct MstConfig
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int P = P_;
};

using MstConfigs
  = ::testing::Types<MstConfig<T8_ECLASS_LINE, 2>, MstConfig<T8_ECLASS_LINE, 3>, MstConfig<T8_ECLASS_QUAD, 2>,
                     MstConfig<T8_ECLASS_QUAD, 3>, MstConfig<T8_ECLASS_HEX, 2>, MstConfig<T8_ECLASS_TRIANGLE, 1>,
                     MstConfig<T8_ECLASS_TRIANGLE, 2>, MstConfig<T8_ECLASS_TRIANGLE, 3>,
                     MstConfig<T8_ECLASS_TRIANGLE, 4>>;

template <typename Config>
class mra_mst: public ::testing::Test {
 public:
  static constexpr t8_eclass Shape = Config::Shape;
  static constexpr unsigned short P = Config::P;

  using element_t = t8_mra::element_data<Shape, 2, P>;
  using mst_t = t8_mra::mst<element_t>;
  using detail_t = typename mst_t::detail_t;
  using lmi_t = t8_mra::levelmultiindex<Shape>;

  static constexpr unsigned int NUM_CHILDREN = lmi_t::NUM_CHILDREN;
  static constexpr unsigned int U = element_t::U_DIM;
  static constexpr unsigned int DOF = element_t::DOF;

  std::vector<t8_mra::mat> mask;

  void
  SetUp () override
  {
    t8_mra::compute_mask<Shape, P> (mask);
  }

  static element_t
  make_leaf (std::size_t seed, double vol = 1.0)
  {
    element_t e;
    e.vol = vol;
    e.order = { 0, 1, 2 };
    for (unsigned int u = 0; u < U; ++u)
      for (unsigned int i = 0; i < DOF; ++i)
        e.u_coeffs[element_t::dg_idx (u, i)] = std::sin (0.7 * (seed + 1) + 1.3 * u + 2.1 * i);
    return e;
  }
};
TYPED_TEST_SUITE (mra_mst, MstConfigs);

/* decomposition then inverse returns the original single family. */
TYPED_TEST (mra_mst, round_trip_identity_one_level)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;

  t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
  t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);

  const auto base = lmi_t (0);
  const auto kids = lmi_t::children (base);
  std::array<element_t, NUM_CHILDREN> orig;
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k) {
    orig[k] = TestFixture::make_leaf (k);
    lmi_map.insert (kids[k], orig[k]);
  }

  mst_t::multiscale_decomposition (0, 1, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (1), 0u) << "family not collapsed";
  ASSERT_TRUE (lmi_map.contains (base));

  mst_t::inverse_multiscale_transformation (0, 1, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (1), NUM_CHILDREN);

  for (unsigned int k = 0; k < NUM_CHILDREN; ++k) {
    const auto &got = lmi_map.get (kids[k]);
    EXPECT_NEAR (got.vol, orig[k].vol, eps);
    for (unsigned int u = 0; u < U; ++u)
      for (unsigned int i = 0; i < DOF; ++i)
        EXPECT_NEAR (got.u_coeffs[element_t::dg_idx (u, i)], orig[k].u_coeffs[element_t::dg_idx (u, i)], eps)
          << "child " << k << " comp " << u << " dof " << i;
  }
}

/* Two-level cascade: collapse NC^2 leaves to the base, reconstruct, compare. */
TYPED_TEST (mra_mst, round_trip_identity_two_levels)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;

  t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
  t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);

  const auto base = lmi_t (0);
  const auto kids = lmi_t::children (base);

  std::vector<lmi_t> leaves;
  std::vector<element_t> orig;
  std::size_t seed = 0;
  for (unsigned int a = 0; a < NUM_CHILDREN; ++a)
    for (const auto &g : lmi_t::children (kids[a])) {
      const auto e = TestFixture::make_leaf (seed++);
      lmi_map.insert (g, e);
      leaves.push_back (g);
      orig.push_back (e);
    }

  mst_t::multiscale_decomposition (0, 2, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (2), 0u);
  ASSERT_EQ (lmi_map.size (1), 0u);
  ASSERT_TRUE (lmi_map.contains (base));

  mst_t::inverse_multiscale_transformation (0, 2, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (2), leaves.size ());

  for (std::size_t n = 0; n < leaves.size (); ++n) {
    const auto &got = lmi_map.get (leaves[n]);
    for (unsigned int u = 0; u < U; ++u)
      for (unsigned int i = 0; i < DOF; ++i)
        EXPECT_NEAR (got.u_coeffs[element_t::dg_idx (u, i)], orig[n].u_coeffs[element_t::dg_idx (u, i)], eps)
          << "leaf " << n << " comp " << u << " dof " << i;
  }
}

/* Conservation: the integral of the parent equals the sum of the children
 * integrals. The cell integral is mean * vol; the constant-mode coefficient
 * gives the mean as u_const for the cartesian basis but u_const / sqrt(vol)
 * for the triangle (Dubiner) basis, so the per-cell integral is
 *   Q = u_const * (cartesian ? vol : sqrt(vol)). */
TYPED_TEST (mra_mst, two_scale_conserves_mass)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  constexpr auto Shape = TestFixture::Shape;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto U = TestFixture::U;

  const auto integral
    = [] (double u_const, double vol) { return u_const * (t8_mra::is_cartesian<Shape> ? vol : std::sqrt (vol)); };

  std::array<element_t, NUM_CHILDREN> sib;
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
    sib[k] = TestFixture::make_leaf (k);

  detail_t parent;
  mst_t::two_scale_family (sib, parent, this->mask);

  for (unsigned int u = 0; u < U; ++u) {
    double child_int = 0.0;
    for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
      child_int += integral (sib[k].u_coeffs[element_t::dg_idx (u, 0)], sib[k].vol);
    const double parent_int = integral (parent.u_coeffs[element_t::dg_idx (u, 0)], parent.vol);
    EXPECT_NEAR (parent_int, child_int, eps) << "component " << u;
  }
}

/* Wavelet kernel: data that lives exactly on the coarse space (a parent with
 * zero details, prolonged to the children) produces zero details when
 * re-analysed, and the parent is recovered. */
TYPED_TEST (mra_mst, details_vanish_on_coarse_representable_data)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;

  t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
  t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);

  const auto base = lmi_t (0);

  /// A coarse parent with arbitrary scaling coefficients and no detail.
  detail_t parent;
  parent.vol = NUM_CHILDREN;  // children get vol 1
  parent.order = { 0, 1, 2 };
  for (unsigned int u = 0; u < U; ++u)
    for (unsigned int i = 0; i < DOF; ++i)
      parent.u_coeffs[element_t::dg_idx (u, i)] = std::cos (0.4 + 1.1 * u + 0.9 * i);

  lmi_map.insert (base, static_cast<const element_t &> (parent));
  d_map.insert (base, parent);  // d_coeffs default to zero

  /// Prolong to the children: they now represent a pure coarse function.
  mst_t::inverse_multiscale_transformation (0, 1, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (1), NUM_CHILDREN);

  /// Re-analyse: details must vanish, parent scaling coefficients recovered.
  mst_t::multiscale_decomposition (0, 1, &lmi_map, d_map, this->mask);
  ASSERT_TRUE (lmi_map.contains (base));

  const auto &re = d_map.get (base);
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
    for (unsigned int u = 0; u < U; ++u)
      for (unsigned int i = 0; i < DOF; ++i)
        EXPECT_NEAR (re.d_coeffs[detail_t::wavelet_idx (k, u, i)], 0.0, eps)
          << "spurious detail at child " << k << " comp " << u << " dof " << i;

  const auto &got = lmi_map.get (base);
  for (unsigned int u = 0; u < U; ++u)
    for (unsigned int i = 0; i < DOF; ++i)
      EXPECT_NEAR (got.u_coeffs[element_t::dg_idx (u, i)], parent.u_coeffs[element_t::dg_idx (u, i)], eps);
}

/* Vanishing moments */
TYPED_TEST (mra_mst, details_vanish_for_projected_polynomial)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  constexpr auto Shape = TestFixture::Shape;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto P = TestFixture::P;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;
  constexpr unsigned int DIM = element_t::DIM;

  using basis_t = t8_mra::basis<Shape, P>;
  const t8_mra::quadrature<Shape> quad (t8_mra::mask_quad_param<Shape, P>);
  const auto children = t8_mra::child_maps<Shape> ();
  constexpr double ref_volume = t8_mra::is_cartesian<Shape> ? 1.0 : 0.5;

  /// Affine form on the parent reference element, raised to P-1: total degree
  /// P-1.
  const auto poly = [] (const std::array<double, DIM> &eta, unsigned int u) {
    const double coeff[3] = { 0.6, 0.4, 0.2 };
    double a = 0.3;
    for (unsigned int d = 0; d < DIM; ++d)
      a += coeff[d] * eta[d];
    return static_cast<double> (u + 1) * std::pow (a, static_cast<int> (P) - 1);
  };

  /// Project the global polynomial onto each child's basis
  std::array<element_t, NUM_CHILDREN> sib;
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k) {
    sib[k].vol = 1.0;
    sib[k].order = { 0, 1, 2 };
    for (unsigned int u = 0; u < U; ++u) {
      std::array<double, DOF> c {};
      for (std::size_t q = 0; q < quad.num_points; ++q) {
        std::array<double, DIM> xi {};
        for (unsigned int d = 0; d < DIM; ++d)
          xi[d] = quad.points[DIM * q + d];
        const auto phi = basis_t::eval (xi);
        const double f = poly (children[k](xi), u);
        const double w = quad.weights[q];
        for (unsigned int i = 0; i < DOF; ++i)
          c[i] += ref_volume * w * f * phi[i];
      }
      for (unsigned int i = 0; i < DOF; ++i)
        sib[k].u_coeffs[element_t::dg_idx (u, i)] = c[i];
    }
  }

  detail_t parent;
  mst_t::two_scale_family (sib, parent, this->mask);

  for (double d : parent.d_coeffs)
    EXPECT_NEAR (d, 0.0, eps) << "non-vanishing detail for a degree <= P-1 polynomial";
}

/* Order-P cancellation */
TYPED_TEST (mra_mst, details_scale_as_h_to_the_P)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  constexpr auto Shape = TestFixture::Shape;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto P = TestFixture::P;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;
  constexpr unsigned int DIM = element_t::DIM;

  using basis_t = t8_mra::basis<Shape, P>;
  const t8_mra::quadrature<Shape> quad (t8_mra::mask_quad_param<Shape, P>);
  const auto children = t8_mra::child_maps<Shape> ();
  constexpr double ref_volume = t8_mra::is_cartesian<Shape> ? 1.0 : 0.5;

  /// f(x) = (sum_d x_d)^P : degree exactly P.
  const auto f = [] (const std::array<double, DIM> &x) {
    double s = 0.0;
    for (unsigned int d = 0; d < DIM; ++d)
      s += x[d];
    return std::pow (s, static_cast<int> (P));
  };

  /// L2 detail norm of f projected onto a family on a parent of physical size h.
  const auto detail_norm = [&] (double h) {
    std::array<element_t, NUM_CHILDREN> sib;
    for (unsigned int k = 0; k < NUM_CHILDREN; ++k) {
      sib[k].vol = 1.0;
      sib[k].order = { 0, 1, 2 };
      for (unsigned int u = 0; u < U; ++u) {
        std::array<double, DOF> c {};
        for (std::size_t q = 0; q < quad.num_points; ++q) {
          std::array<double, DIM> xi {};
          for (unsigned int d = 0; d < DIM; ++d)
            xi[d] = quad.points[DIM * q + d];
          const auto phi = basis_t::eval (xi);
          const auto eta = children[k](xi);  // child ref -> parent ref
          std::array<double, DIM> x {};
          for (unsigned int d = 0; d < DIM; ++d)
            x[d] = h * eta[d];  // parent ref -> physical (size h)
          const double val = f (x);
          const double w = quad.weights[q];
          for (unsigned int i = 0; i < DOF; ++i)
            c[i] += ref_volume * w * val * phi[i];
        }
        for (unsigned int i = 0; i < DOF; ++i)
          sib[k].u_coeffs[element_t::dg_idx (u, i)] = c[i];
      }
    }
    detail_t parent;
    mst_t::two_scale_family (sib, parent, this->mask);
    double s = 0.0;
    for (double d : parent.d_coeffs)
      s += d * d;
    return std::sqrt (s);
  };

  const double D1 = detail_norm (1.0);
  const double D2 = detail_norm (0.5);
  EXPECT_GT (D1, eps) << "a degree-P polynomial must leave a detail (it is not in the coarse space)";
  EXPECT_NEAR (D2 / D1, std::pow (0.5, static_cast<int> (P)), eps)
    << "details must decay as h^P (order-P cancellation)";
}

/* Round trip on a graded (adaptive) grid: one branch refined a level deeper
 * than the rest. */
TYPED_TEST (mra_mst, round_trip_identity_graded_grid)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;

  t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
  t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);

  const auto base = lmi_t (0);
  const auto kids = lmi_t::children (base);  // level 1

  std::vector<lmi_t> leaves;
  std::vector<element_t> orig;
  std::size_t seed = 0;

  /// kids[0] is refined one level deeper: its NUM_CHILDREN children are level-2 leaves.
  for (const auto &g : lmi_t::children (kids[0])) {
    const auto e = TestFixture::make_leaf (seed++, 1.0);
    lmi_map.insert (g, e);
    leaves.push_back (g);
    orig.push_back (e);
  }

  /// The remaining children of the base stay as level-1 leaves (vol NUM_CHILDREN times the
  /// level-2 cells, the consistent grading).
  for (unsigned int k = 1; k < NUM_CHILDREN; ++k) {
    const auto e = TestFixture::make_leaf (seed++, static_cast<double> (NUM_CHILDREN));
    lmi_map.insert (kids[k], e);
    leaves.push_back (kids[k]);
    orig.push_back (e);
  }

  mst_t::multiscale_decomposition (0, 2, &lmi_map, d_map, this->mask);
  ASSERT_TRUE (lmi_map.contains (base));
  ASSERT_EQ (lmi_map.size (), 1u) << "graded grid not fully collapsed to the base";

  mst_t::inverse_multiscale_transformation (0, 2, &lmi_map, d_map, this->mask);
  ASSERT_EQ (lmi_map.size (), leaves.size ());

  for (std::size_t n = 0; n < leaves.size (); ++n) {
    ASSERT_TRUE (lmi_map.contains (leaves[n])) << "missing leaf " << n;
    const auto &got = lmi_map.get (leaves[n]);
    for (unsigned int u = 0; u < U; ++u)
      for (unsigned int i = 0; i < DOF; ++i)
        EXPECT_NEAR (got.u_coeffs[element_t::dg_idx (u, i)], orig[n].u_coeffs[element_t::dg_idx (u, i)], eps)
          << "leaf " << n << " comp " << u << " dof " << i;
  }
}

/* Sanity guard for the vanishing test: generic (non-coarse) data must produce
 * at least one non-trivial detail */
TYPED_TEST (mra_mst, generic_data_produces_nonzero_details)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;

  if constexpr (TestFixture::DOF == 1) {
    GTEST_SKIP () << "P=1 has only the constant mode; no detail can be non-zero";
  }
  else {
    std::array<element_t, NUM_CHILDREN> sib;
    for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
      sib[k] = TestFixture::make_leaf (k);

    detail_t parent;
    mst_t::two_scale_family (sib, parent, this->mask);

    double max_detail = 0.0;
    for (double d : parent.d_coeffs)
      max_detail = std::max (max_detail, std::abs (d));
    EXPECT_GT (max_detail, eps);
  }
}

/* The non-destructive analysis keeps the leaves in lmi_map and writes the
 * parent (scaling coeffs + details) into d_map; the details equal a direct
 * two_scale_family of the same children. */
TYPED_TEST (mra_mst, multiscale_transformation_is_non_destructive)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  constexpr auto DOF = TestFixture::DOF;
  constexpr auto U = TestFixture::U;

  t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
  t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);

  const auto base = lmi_t (0);
  const auto kids = lmi_t::children (base);
  std::array<element_t, NUM_CHILDREN> sib;
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k) {
    sib[k] = TestFixture::make_leaf (k);
    lmi_map.insert (kids[k], sib[k]);
  }

  mst_t::multiscale_transformation (0, 1, &lmi_map, d_map, this->mask);

  /// Leaves untouched, parent not inserted.
  EXPECT_EQ (lmi_map.size (1), NUM_CHILDREN);
  EXPECT_FALSE (lmi_map.contains (base));
  for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
    EXPECT_TRUE (lmi_map.contains (kids[k]));

  /// d_map holds the parent, equal to a direct two-scale of the same family.
  ASSERT_TRUE (d_map.contains (base));
  detail_t direct;
  mst_t::two_scale_family (sib, direct, this->mask);
  const auto &got = d_map.get (base);
  for (unsigned int u = 0; u < U; ++u)
    for (unsigned int i = 0; i < DOF; ++i) {
      EXPECT_NEAR (got.u_coeffs[element_t::dg_idx (u, i)], direct.u_coeffs[element_t::dg_idx (u, i)], eps);
      for (unsigned int k = 0; k < NUM_CHILDREN; ++k)
        EXPECT_NEAR (got.d_coeffs[detail_t::wavelet_idx (k, u, i)], direct.d_coeffs[detail_t::wavelet_idx (k, u, i)],
                     eps);
    }
}

/* An incomplete family (a missing sibling) is skipped by both the analysis and
 * the destructive collapse: no parent appears, the present leaves stay. */
TYPED_TEST (mra_mst, incomplete_family_is_skipped)
{
  using element_t = typename TestFixture::element_t;
  using detail_t = typename TestFixture::detail_t;
  using mst_t = typename TestFixture::mst_t;
  using lmi_t = typename TestFixture::lmi_t;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;

  const auto base = lmi_t (0);
  const auto kids = lmi_t::children (base);

  // Only NUM_CHILDREN-1 of the NUM_CHILDREN siblings are present.
  auto build = [&] (t8_mra::levelindex_map<lmi_t, element_t> &lmi_map) {
    for (unsigned int k = 0; k < NUM_CHILDREN - 1; ++k)
      lmi_map.insert (kids[k], TestFixture::make_leaf (k));
  };

  {
    t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
    t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);
    build (lmi_map);
    mst_t::multiscale_transformation (0, 1, &lmi_map, d_map, this->mask);
    EXPECT_FALSE (d_map.contains (base)) << "incomplete family produced a detail";
    EXPECT_EQ (lmi_map.size (1), NUM_CHILDREN - 1);
  }
  {
    t8_mra::levelindex_map<lmi_t, element_t> lmi_map (3);
    t8_mra::levelindex_map<lmi_t, detail_t> d_map (3);
    build (lmi_map);
    mst_t::multiscale_decomposition (0, 1, &lmi_map, d_map, this->mask);
    EXPECT_FALSE (lmi_map.contains (base)) << "incomplete family collapsed";
    EXPECT_EQ (lmi_map.size (1), NUM_CHILDREN - 1) << "present leaves must remain";
  }
}

/* Scaling policy values per shape. */
TYPED_TEST (mra_mst, scaling_policy_values)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr auto NUM_CHILDREN = TestFixture::NUM_CHILDREN;
  using policy = t8_mra::mst_scaling_policy<Shape>;

  EXPECT_EQ (policy::inverse_scaling_factor (), 1.0);
  if constexpr (t8_mra::is_cartesian<Shape>)
    EXPECT_EQ (policy::forward_scaling_factor (NUM_CHILDREN), 1.0 / static_cast<double> (NUM_CHILDREN));
  else
    EXPECT_EQ (policy::forward_scaling_factor (NUM_CHILDREN), 1.0);
}

}  // namespace

#endif  // T8_ENABLE_MRA

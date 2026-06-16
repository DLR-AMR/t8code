#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/core/shape_traits.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

#include <array>
#include <cmath>
#include <vector>

namespace
{

template <t8_eclass TShape, int P_>
struct NumConfig
{
  static constexpr t8_eclass Shape = TShape;
  static constexpr int P = P_;
  static constexpr int DIM = t8_mra::shape_traits<TShape>::DIM;
};

using NumConfigs
  = ::testing::Types<NumConfig<T8_ECLASS_LINE, 2>, NumConfig<T8_ECLASS_LINE, 4>, NumConfig<T8_ECLASS_QUAD, 2>,
                     NumConfig<T8_ECLASS_QUAD, 3>, NumConfig<T8_ECLASS_HEX, 2>, NumConfig<T8_ECLASS_TRIANGLE, 1>,
                     NumConfig<T8_ECLASS_TRIANGLE, 2>, NumConfig<T8_ECLASS_TRIANGLE, 3>,
                     NumConfig<T8_ECLASS_TRIANGLE, 4>>;

/// Volume of the reference cell: 1 for the cartesian [0,1]^DIM, 1/2 for the
/// reference triangle. The quadrature weights are normalized to sum to 1 on
/// every shape, so the integral of f over the reference cell is
/// ref_volume * sum_q w_q f(x_q).
template <t8_eclass Shape>
inline constexpr double ref_volume = t8_mra::is_cartesian<Shape> ? 1.0 : 0.5;

/// A handful of points strictly inside the reference cell (away from edges so
/// finite differences stay valid).
template <t8_eclass Shape, int DIM>
std::vector<std::array<double, DIM>>
interior_points ()
{
  std::vector<std::array<double, DIM>> pts;
  if constexpr (t8_mra::is_cartesian<Shape>) {
    std::array<double, DIM> a;
    a.fill (0.3);
    std::array<double, DIM> b;
    b.fill (0.55);
    pts.push_back (a);
    pts.push_back (b);
    if constexpr (DIM >= 2)
      pts.push_back ([] {
        std::array<double, DIM> c;
        c.fill (0.4);  // all dims interior; FD steps must stay inside [0,1]
        c[0] = 0.2;
        c[1] = 0.7;
        return c;
      }());
  }
  else {
    pts.push_back ({ 0.25, 0.25 });
    pts.push_back ({ 0.5, 0.3 });
    pts.push_back ({ 0.2, 0.6 });
  }
  return pts;
}

template <typename Config>
class mra_num: public ::testing::Test {
 public:
  static constexpr t8_eclass Shape = Config::Shape;
  static constexpr int P = Config::P;
  static constexpr int DIM = Config::DIM;
  static constexpr int DOF = t8_mra::shape_traits<Shape>::dof (P);

  using basis_t = t8_mra::basis<Shape, P>;

  static constexpr double eps = 1e-12;

  // Exact to degree >= 2(P-1) on every shape, enough for the mass matrix.
  t8_mra::quadrature<Shape> quad { t8_mra::mask_quad_param<Shape, P> };
};
TYPED_TEST_SUITE (mra_num, NumConfigs);

/* Weights are positive, sum to one (normalized rule) */
TYPED_TEST (mra_num, quadrature_weights_and_points)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int DIM = TestFixture::DIM;
  constexpr auto eps = TestFixture::eps;
  const auto &quad = this->quad;

  double wsum = 0.0;
  for (std::size_t i = 0; i < quad.num_points; ++i) {
    EXPECT_GT (quad.weights[i], 0.0);
    wsum += quad.weights[i];

    if constexpr (t8_mra::is_cartesian<Shape>) {
      for (int d = 0; d < DIM; ++d) {
        EXPECT_GE (quad.points[DIM * i + d], 0.0);
        EXPECT_LE (quad.points[DIM * i + d], 1.0);
      }
    }
    else {
      const double x = quad.points[2 * i + 0];
      const double y = quad.points[2 * i + 1];
      EXPECT_GE (x, 0.0);
      EXPECT_GE (y, 0.0);
      EXPECT_LE (x + y, 1.0 + eps);
    }
  }
  EXPECT_NEAR (wsum, 1.0, eps);
}

/* Polynomial exactness: cartesian integrates x0^(2n-1) = 1/(2n) exactly;
 * triangle integrates x^a y^b = a! b! / (a+b+2)! exactly up to the rule degree. */
TYPED_TEST (mra_num, quadrature_is_exact)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int DIM = TestFixture::DIM;
  constexpr int P = TestFixture::P;
  constexpr auto eps = TestFixture::eps;
  const auto &quad = this->quad;

  if constexpr (t8_mra::is_cartesian<Shape>) {
    const int n = t8_mra::mask_quad_param<Shape, P>;  // points per axis
    const int m = 2 * n - 1;                          // highest exactly integrable degree

    // Single-axis monomial x0^m.
    double axis = 0.0;
    for (std::size_t i = 0; i < quad.num_points; ++i)
      axis += quad.weights[i] * std::pow (quad.points[DIM * i + 0], m);
    EXPECT_NEAR (axis, 1.0 / (m + 1), eps);

    // Full-dimensional monomial prod_d x_d^m: exercises the tensor weight
    // product across every axis. Exact value (1/(m+1))^DIM.
    double mixed = 0.0;
    for (std::size_t i = 0; i < quad.num_points; ++i) {
      double f = quad.weights[i];
      for (int d = 0; d < DIM; ++d)
        f *= std::pow (quad.points[DIM * i + d], m);
      mixed += f;
    }
    EXPECT_NEAR (mixed, std::pow (1.0 / (m + 1), DIM), eps);
  }
  else {
    const auto mono = [&] (int a, int b) {
      double s = 0.0;
      for (std::size_t i = 0; i < quad.num_points; ++i)
        s += quad.weights[i] * std::pow (quad.points[2 * i + 0], a) * std::pow (quad.points[2 * i + 1], b);
      return s;
    };
    const auto fact = [] (int k) {
      double f = 1.0;
      for (int i = 2; i <= k; ++i)
        f *= i;
      return f;
    };
    const auto exact = [&] (int a, int b) { return fact (a) * fact (b) / fact (a + b + 2); };

    // integral over the reference triangle = ref_volume * sum_q w_q f(x_q)
    EXPECT_NEAR (ref_volume<Shape> * mono (0, 0), exact (0, 0), eps);  // area = 1/2
    EXPECT_NEAR (ref_volume<Shape> * mono (1, 0), exact (1, 0), eps);  // 1/6
    EXPECT_NEAR (ref_volume<Shape> * mono (1, 1), exact (1, 1), eps);  // 1/24
  }
}

/* The reference basis is orthonormal */
TYPED_TEST (mra_num, basis_is_orthonormal)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int DIM = TestFixture::DIM;
  constexpr int DOF = TestFixture::DOF;
  constexpr auto eps = TestFixture::eps;
  using basis_t = typename TestFixture::basis_t;
  const auto &quad = this->quad;

  for (int i = 0; i < DOF; ++i)
    for (int j = 0; j < DOF; ++j) {
      double mij = 0.0;

      for (std::size_t k = 0; k < quad.num_points; ++k) {
        std::array<double, DIM> x {};
        for (int d = 0; d < DIM; ++d)
          x[d] = quad.points[DIM * k + d];
        const auto phi = basis_t::eval (x);
        mij += quad.weights[k] * phi[i] * phi[j];
      }

      mij *= ref_volume<Shape>;  // integral over the reference cell
      EXPECT_NEAR (mij, (i == j) ? 1.0 : 0.0, eps) << "mass(" << i << "," << j << ")";
    }
}

}  // namespace

#endif  // T8_ENABLE_MRA

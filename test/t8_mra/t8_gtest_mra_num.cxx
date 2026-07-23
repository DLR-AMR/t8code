#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

#include <t8_mra/core/shape_traits.hxx>
#include <t8_mra/num/basis/basis.hxx>
#include <t8_mra/num/geometry.hxx>
#include <t8_mra/num/mask_coefficients.hxx>
#include <t8_mra/num/mat.hxx>
#include <t8_mra/num/nodal_to_modal.hxx>
#include <t8_mra/num/quadrature/quadrature.hxx>

#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace
{

constexpr double eps = 1e-12;

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

/// Volume of the reference cell:
///   1 for the cartesian [0,1]^DIM
///   1/2 for the reference triangle
template <t8_eclass Shape>
inline constexpr double ref_volume = t8_mra::is_cartesian<Shape> ? 1.0 : 0.5;

/// Points strictly inside the reference cell
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
        c.fill (0.4);
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

  /// Exact to degree >= 2(P-1) on every shape
  t8_mra::quadrature<Shape> quad { t8_mra::mask_quad_param<Shape, P> };
};
TYPED_TEST_SUITE (mra_num, NumConfigs);

/* Weights are positive, sum to one (normalized rule) */
TYPED_TEST (mra_num, quadrature_weights_and_points)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int DIM = TestFixture::DIM;
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
  const auto &quad = this->quad;

  if constexpr (t8_mra::is_cartesian<Shape>) {
    const int n = t8_mra::mask_quad_param<Shape, P>;  // points per axis
    const int m = 2 * n - 1;                          // highest exactly integrable degree

    /// Single-axis monomial x0^m.
    double axis = 0.0;
    for (std::size_t i = 0; i < quad.num_points; ++i)
      axis += quad.weights[i] * std::pow (quad.points[DIM * i + 0], m);
    EXPECT_NEAR (axis, 1.0 / (m + 1), eps);

    /// Full-dimensional monomial prod_d x_d^m: exercises the tensor weight
    /// product across every axis. Exact value (1/(m+1))^DIM.
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

    /// integral over the reference triangle = ref_volume * sum_q w_q f(x_q)
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

/* eval_gradient matches the analytic gradient of a known polynomial.
 *
 * The basis spans every polynomial of degree <= P-1, so the monomial sum
 *   u(x) = sum_d x_d^(P-1)
 * is represented exactly by its projection coefficients c_i = integral u phi_i.
 */
TYPED_TEST (mra_num, gradient_matches_exact_polynomial)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int DIM = TestFixture::DIM;
  constexpr int DOF = TestFixture::DOF;
  constexpr int P = TestFixture::P;
  using basis_t = typename TestFixture::basis_t;
  const auto &quad = this->quad;

  const auto u = [] (const std::array<double, DIM> &x) {
    double s = 0.0;
    for (int d = 0; d < DIM; ++d)
      s += std::pow (x[d], P - 1);
    return s;
  };
  const auto grad_u
    = [] (const std::array<double, DIM> &x, int dir) { return (P >= 2) ? (P - 1) * std::pow (x[dir], P - 2) : 0.0; };

  std::array<double, DOF> c {};
  for (std::size_t q = 0; q < quad.num_points; ++q) {
    std::array<double, DIM> xq {};

    for (int d = 0; d < DIM; ++d)
      xq[d] = quad.points[DIM * q + d];

    const auto phi = basis_t::eval (xq);
    const double uw = ref_volume<Shape> * quad.weights[q] * u (xq);

    for (int i = 0; i < DOF; ++i)
      c[i] += uw * phi[i];
  }

  for (const auto &x : interior_points<Shape, DIM> ()) {
    const auto grad = basis_t::eval_gradient (x);
    for (int dir = 0; dir < DIM; ++dir) {
      double recon = 0.0;
      for (int i = 0; i < DOF; ++i)
        recon += c[i] * grad[dir][i];
      EXPECT_NEAR (recon, grad_u (x, dir), eps) << "dir " << dir;
    }
  }
}

/* Two-scale refinement equation: the parent basis evaluated at the child-mapped
 * point equals a mask-weighted sum of child basis functions,
 *   phi_j(Phi_k xi) = (ref_vol / mask_norm) * sum_i mask_k(i,j) phi_i(xi).
 */
TYPED_TEST (mra_num, mask_satisfies_refinement_equation)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int P = TestFixture::P;
  constexpr int DIM = TestFixture::DIM;
  constexpr int DOF = TestFixture::DOF;
  using basis_t = typename TestFixture::basis_t;

  std::vector<t8_mra::mat> mask;
  t8_mra::compute_mask<Shape, P> (mask);

  const auto children = t8_mra::child_maps<Shape> ();
  const double factor = ref_volume<Shape> / t8_mra::mask_norm<Shape>;

  for (std::size_t k = 0; k < children.size (); ++k)
    for (const auto &xi : interior_points<Shape, DIM> ()) {
      const auto phi_child = basis_t::eval (xi);
      const auto phi_parent = basis_t::eval (children[k](xi));
      for (int j = 0; j < DOF; ++j) {
        double rhs = 0.0;
        for (int i = 0; i < DOF; ++i)
          rhs += mask[k](i, j) * phi_child[i];
        EXPECT_NEAR (phi_parent[j], factor * rhs, eps) << "child " << k << " parent dof " << j;
      }
    }
}

/* ---- geometry: reference -> physical cartesian map ---- */
TEST (mra_geometry, deref_1d_maps_endpoints_and_midpoint)
{
  EXPECT_NEAR (t8_mra::deref_1d (0.0, 2.0, 5.0), 2.0, eps);
  EXPECT_NEAR (t8_mra::deref_1d (1.0, 2.0, 5.0), 5.0, eps);
  EXPECT_NEAR (t8_mra::deref_1d (0.5, 2.0, 5.0), 3.5, eps);
}

TEST (mra_geometry, deref_maps_box_corners)
{
  const std::array<double, 2> lo { -1.0, 3.0 }, hi { 1.0, 4.0 };
  EXPECT_EQ (t8_mra::deref<2> ({ 0.0, 0.0 }, lo, hi), lo);
  EXPECT_EQ (t8_mra::deref<2> ({ 1.0, 1.0 }, lo, hi), hi);
  const auto mid = t8_mra::deref<2> ({ 0.5, 0.5 }, lo, hi);
  EXPECT_NEAR (mid[0], 0.0, eps);
  EXPECT_NEAR (mid[1], 3.5, eps);
}

/* extract_cartesian_vertices picks the lower corner (index 0) and the upper
 * corner (index 2 for a 2D cell) out of the t8code vertex layout. */
TEST (mra_geometry, extract_cartesian_vertices_picks_min_max)
{
  /// QUAD physical vertices
  const double verts[4][3] = { { 1.0, 2.0, 0.0 }, { 3.0, 2.0, 0.0 }, { 4.0, 5.0, 0.0 }, { 3.0, 5.0, 0.0 } };
  std::array<double, 2> lo {}, hi {};

  t8_mra::extract_cartesian_vertices<2> (verts, lo, hi);
  EXPECT_EQ (lo, (std::array<double, 2> { 1.0, 2.0 }));
  EXPECT_EQ (hi, (std::array<double, 2> { 4.0, 5.0 }));
}

/* transform_quad_points maps flattened reference points onto the physical box. */
TEST (mra_geometry, transform_quad_points_maps_onto_box)
{
  const std::array<double, 2> lo { 1.0, 2.0 }, hi { 4.0, 5.0 };
  const std::vector<double> ref { 0.0, 0.0, 1.0, 1.0, 0.5, 0.5 };  // 3 points
  const auto phys = t8_mra::transform_quad_points<2> (ref, 3, lo, hi);

  ASSERT_EQ (phys.size (), 6u);
  EXPECT_NEAR (phys[0], 1.0, eps);  // (0,0) -> lo
  EXPECT_NEAR (phys[1], 2.0, eps);
  EXPECT_NEAR (phys[2], 4.0, eps);  // (1,1) -> hi
  EXPECT_NEAR (phys[3], 5.0, eps);
  EXPECT_NEAR (phys[4], 2.5, eps);  // (0.5,0.5) -> centre
  EXPECT_NEAR (phys[5], 3.5, eps);
}

/* ---- mat: dense matrix and its LU solver ---- */
TEST (mra_mat, element_access_and_resize)
{
  t8_mra::mat A (2, 3);
  EXPECT_EQ (A.rows (), 2u);
  EXPECT_EQ (A.cols (), 3u);

  A (1, 2) = 7.0;
  A (0, 0) = -1.0;
  EXPECT_EQ (A (1, 2), 7.0);
  EXPECT_EQ (A (0, 0), -1.0);
  EXPECT_EQ (A (0, 1), 0.0);  // value-initialized
#if T8_ENABLE_DEBUG
  EXPECT_THROW (A (2, 0), std::out_of_range);  // bounds check is debug-only
#endif

  A.resize (4, 4);
  EXPECT_EQ (A.rows (), 4u);
  EXPECT_EQ (A.cols (), 4u);
  EXPECT_EQ (A (3, 3), 0.0);  // resize clears
}

/* LU factor + solve recovers the known solution of A x = b (with pivoting). */
TEST (mra_mat, lu_solve_recovers_known_solution)
{
  /// Non-symmetric, well-conditioned 3x3; row order forces a pivot.
  t8_mra::mat A (3, 3, { 0.0, 2.0, 1.0, 1.0, 3.0, -1.0, 2.0, 1.0, 4.0 });
  const std::vector<double> x_true { 1.0, -2.0, 3.0 };

  std::vector<double> b (3, 0.0);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      b[i] += A (i, j) * x_true[j];

  std::vector<size_t> p;
  std::vector<double> x = b;
  t8_mra::lu_factors (A, p);  // A becomes its LU factors, p the pivot order
  t8_mra::lu_solve (A, p, x);

  for (size_t i = 0; i < 3; ++i)
    EXPECT_NEAR (x[i], x_true[i], eps) << "component " << i;
}

/* ---- nodal <-> modal converters ---- */
using ConvConfigs
  = ::testing::Types<NumConfig<T8_ECLASS_LINE, 2>, NumConfig<T8_ECLASS_LINE, 3>, NumConfig<T8_ECLASS_QUAD, 2>,
                     NumConfig<T8_ECLASS_QUAD, 3>, NumConfig<T8_ECLASS_HEX, 2>>;

template <typename Config>
class mra_nodal_modal: public ::testing::Test {
 public:
  static constexpr t8_eclass Shape = Config::Shape;
  static constexpr int P = Config::P;
  static constexpr int DIM = Config::DIM;
  static constexpr int DOF = t8_mra::shape_traits<Shape>::dof (P);

  using basis_t = t8_mra::basis<Shape, P>;
  using node_set = std::array<std::array<double, DIM>, DOF>;

  /// Tensor P distinct 1D nodes onto DOF cell nodes (basis-index decomposition,
  /// first coordinate fastest); any distinct set makes the Vandermonde
  /// nonsingular.
  static node_set
  tensor_nodes (const std::array<double, P> &x1d)
  {
    node_set nd {};
    for (int p = 0; p < DOF; ++p) {
      int idx = p;
      for (int d = 0; d < DIM; ++d) {
        nd[p][d] = x1d[idx % P];
        idx /= P;
      }
    }
    return nd;
  }

  /// Uniform (equispaced) and nonuniform (Chebyshev-Lobatto, clustered at the
  /// ends) node grids: the converters must handle both.
  std::vector<node_set> node_sets = {
    tensor_nodes ([] {
      std::array<double, P> x {};
      for (int i = 0; i < P; ++i)
        x[i] = i / static_cast<double> (P - 1);
      return x;
    }()),
    tensor_nodes ([] {
      std::array<double, P> x {};
      for (int i = 0; i < P; ++i)
        x[i] = 0.5 * (1.0 - std::cos (M_PI * i / (P - 1)));
      return x;
    }()),
  };
};
TYPED_TEST_SUITE (mra_nodal_modal, ConvConfigs);

/* modal -> nodal -> modal is the identity (square Vandermonde interpolation),
 * with two components to cover the component-major layout. */
TYPED_TEST (mra_nodal_modal, roundtrip_recovers_modal_coefficients)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int P = TestFixture::P;
  constexpr int DOF = TestFixture::DOF;
  constexpr unsigned int U = 2;

  for (const auto &nodes : this->node_sets) {
    const t8_mra::nodal_to_modal<Shape, U, P> to_modal (nodes);
    const t8_mra::modal_to_nodal<Shape, U, P> to_nodal (nodes);

    std::array<double, U * DOF> modal {};
    for (unsigned int i = 0; i < U * DOF; ++i)
      modal[i] = std::sin (0.7 * i + 0.3) + 0.5;

    const auto nodal = to_nodal (std::span<const double> (modal.data (), modal.size ()));
    const auto back = to_modal (std::span<const double> (nodal.data (), nodal.size ()));

    for (unsigned int i = 0; i < U * DOF; ++i)
      EXPECT_NEAR (back[i], modal[i], eps) << "coeff " << i;
  }
}

/* Nodal interpolation is exact for polynomials in the basis span: sampling
 * u(x) = sum_d x_d^(P-1) at the nodes and converting to modal reconstructs u
 * everywhere. */
TYPED_TEST (mra_nodal_modal, interpolation_is_exact_for_polynomials)
{
  constexpr auto Shape = TestFixture::Shape;
  constexpr int P = TestFixture::P;
  constexpr int DIM = TestFixture::DIM;
  constexpr int DOF = TestFixture::DOF;
  using basis_t = typename TestFixture::basis_t;

  const auto u = [] (const std::array<double, DIM> &x) {
    double s = 0.0;
    for (int d = 0; d < DIM; ++d)
      s += std::pow (x[d], P - 1);
    return s;
  };

  for (const auto &nodes : this->node_sets) {
    std::array<double, DOF> nodal {};
    for (int j = 0; j < DOF; ++j)
      nodal[j] = u (nodes[j]);

    const t8_mra::nodal_to_modal<Shape, 1, P> to_modal (nodes);
    const auto modal = to_modal (std::span<const double> (nodal.data (), nodal.size ()));

    for (const auto &x : interior_points<Shape, DIM> ()) {
      const auto phi = basis_t::eval (x);
      double recon = 0.0;
      for (int i = 0; i < DOF; ++i)
        recon += modal[i] * phi[i];  // cartesian normalization is 1
      EXPECT_NEAR (recon, u (x), eps);
    }
  }
}

}  // namespace

#endif  // T8_ENABLE_MRA

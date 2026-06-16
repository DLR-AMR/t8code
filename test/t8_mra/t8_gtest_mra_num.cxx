#include <gtest/gtest.h>

#ifdef T8_ENABLE_MRA

#include <t8.h>
#include <t8_eclass/t8_eclass.h>

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

}  // namespace

#endif  // T8_ENABLE_MRA

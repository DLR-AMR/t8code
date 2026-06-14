#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <concepts>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis_functions.hxx"
#include "t8_mra/num/legendre_basis.hxx"

namespace t8_mra
{

// ============================================================================
// Reference-element polynomial basis
// ============================================================================
// One specialization of basis<TShape, P> per element shape is the single place
// that defines a shape's function space. Everything that evaluates the basis
// (projection in dg_basis, VTK output) goes through this interface.
//
// Adding a new shape: specialize basis<NewShape, P> with DIM/DOF, eval(), and
// normalization() (plus eval_gradient() if its projection needs gradients).

/// Minimal interface every basis specialization provides.
template <typename B>
concept reference_basis = requires (std::array<double, static_cast<std::size_t> (B::DIM)> x, double vol) {
  { B::DIM } -> std::convertible_to<int>;
  { B::DOF } -> std::convertible_to<int>;
  { B::eval (x) } -> std::same_as<std::array<double, static_cast<std::size_t> (B::DOF)>>;
  { B::normalization (vol) } -> std::convertible_to<double>;
};

template <t8_eclass TShape, int P>
struct basis;

/// Cartesian shapes (LINE, QUAD, HEX): tensor product of 1D Legendre modes,
/// the basis index decomposed lexicographically (first coordinate fastest).
/// Orthonormal on the reference cell, so no volume normalization.
template <t8_eclass TShape, int P>
  requires is_cartesian<TShape>
struct basis<TShape, P>
{
  static constexpr int DIM = shape_traits<TShape>::DIM;
  static constexpr int DOF = shape_traits<TShape>::dof (P);

  static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    std::array<double, DOF> res = {};
    for (int p = 0; p < DOF; ++p) {
      double v = 1.0;
      int idx = p;
      for (int d = 0; d < DIM; ++d) {
        v *= phi_1d (x[d], idx % P);
        idx /= P;
      }
      res[p] = v;
    }
    return res;
  }

  static std::array<std::array<double, DOF>, DIM>
  eval_gradient (const std::array<double, DIM> &x)
  {
    std::array<std::array<double, DOF>, DIM> grad = {};
    for (int dir = 0; dir < DIM; ++dir) {
      for (int p = 0; p < DOF; ++p) {
        double v = 1.0;
        int idx = p;
        for (int d = 0; d < DIM; ++d) {
          const int deg = idx % P;
          idx /= P;
          v *= (d == dir) ? phi_prime_1d<P> (x[d], deg) : phi_1d (x[d], deg);
        }
        grad[dir][p] = v;
      }
    }
    return grad;
  }

  static constexpr double
  normalization (double) noexcept
  {
    return 1.0;
  }
};

/// Triangle: the orthonormal Dubiner basis in barycentric coords
/// (x = {lambda0, lambda1}). Orthonormal on the reference triangle (area 1/2),
/// so the physical basis scales by sqrt(1/(2*vol)).
template <int P>
struct basis<T8_ECLASS_TRIANGLE, P>
{
  static constexpr int DIM = 2;
  static constexpr int DOF = shape_traits<T8_ECLASS_TRIANGLE>::dof (P);

  static std::array<double, DOF>
  eval (const std::array<double, DIM> &x)
  {
    std::array<double, DOF> res = {};
    for (int i = 0; i < DOF; ++i)
      res[i] = scaling_function (i, x[0], x[1]);
    return res;
  }

  static double
  normalization (double vol)
  {
    return std::sqrt (1.0 / (2.0 * vol));
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

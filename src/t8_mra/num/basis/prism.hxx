#pragma once

#ifdef T8_ENABLE_MRA

#include <array>

#include <t8_eclass/t8_eclass.h>

#include "t8_mra/core/shape_traits.hxx"
#include "t8_mra/num/basis/dubiner.hxx"
#include "t8_mra/num/basis/legendre.hxx"

namespace t8_mra
{

namespace detail
{
template <int P>
inline constexpr int prism_tri_dof = shape_traits<T8_ECLASS_TRIANGLE>::dof (P);
}  // namespace detail

/// I-th orthonormal prism scaling function: the Dubiner function of the
/// triangle factor times the Legendre mode of the line factor. Orthonormal on
/// the reference prism (volume 1/2). tau1, tau2 barycentric, z in [0,1].
template <int P, int I>
[[nodiscard]] double
prism_scaling_function (double tau1, double tau2, double z)
{
  constexpr int tri_index = I % detail::prism_tri_dof<P>;
  constexpr int z_degree = I / detail::prism_tri_dof<P>;

  return scaling_function<tri_index> (tau1, tau2) * phi_1d (z, z_degree);
}

/// Reference gradient {d/dtau1, d/dtau2, d/dz} of the I-th prism scaling function.
template <int P, int I>
[[nodiscard]] std::array<double, 3>
prism_scaling_function_gradient (double tau1, double tau2, double z)
{
  constexpr int tri_index = I % detail::prism_tri_dof<P>;
  constexpr int z_degree = I / detail::prism_tri_dof<P>;

  const auto tri_grad = scaling_function_gradient<tri_index> (tau1, tau2);
  const double leg = phi_1d (z, z_degree);

  return { tri_grad[0] * leg, tri_grad[1] * leg,
           scaling_function<tri_index> (tau1, tau2) * phi_prime_1d<P> (z, z_degree) };
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

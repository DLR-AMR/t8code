#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>

#include "t8_mra/num/legendre_basis.hxx"

namespace t8_mra
{

/**
 * @brief Multiindex type for D-dimensional indices
 *
 * Used to represent tensor product basis indices.
 * For example, in 2D with P=3: multiindex {1,2} represents phi_1(x) * phi_2(y)
 *
 * @tparam DIM Spatial dimension
 */
template <unsigned int DIM>
using multiindex = std::array<int, DIM>;

/**
 * @brief Computes the number of basis functions for tensor-structured basis
 *
 * For a tensor basis of polynomial degree P-1 in D dimensions,
 * the number of basis functions is P^D.
 *
 * @param P Polynomial order (degree is P-1)
 * @param D Spatial dimension
 * @return constexpr int Number of basis functions
 */
constexpr int
tensor_basis_size (int P, int D)
{
  int result = 1;
  for (int i = 0; i < D; ++i)
    result *= P;

  return result;
}

/**
 * @brief Generates the multiindex set (pset) for tensor-structured basis
 *
 * Creates all possible combinations of 1D basis indices for a tensor product basis.
 * For D dimensions and polynomial order P, generates P^D multiindices.
 *
 * Example for D=2, P=3:
 *   pset[0] = {0,0}, pset[1] = {0,1}, pset[2] = {0,2},
 *   pset[3] = {1,0}, pset[4] = {1,1}, pset[5] = {1,2},
 *   pset[6] = {2,0}, pset[7] = {2,1}, pset[8] = {2,2}
 *
 * The ordering is lexicographic with the last dimension varying fastest.
 *
 * @tparam DIM Spatial dimension
 * @param P Polynomial order (degree is P-1)
 * @return std::vector<multiindex<DIM>> Vector of all multiindices
 */
template <unsigned int DIM>
std::vector<multiindex<DIM>>
generate_tensor_pset (int P)
{
  const int num_basis = tensor_basis_size (P, DIM);
  std::vector<multiindex<DIM>> pset (num_basis);

  if constexpr (DIM == 1) {
    // 1D case: simple enumeration
    for (int i = 0; i < P; ++i) {
      pset[i][0] = i;
    }
  }
  else if constexpr (DIM == 2) {
    int idx = 0;
    for (int j = 0; j < P; ++j) {
      for (int i = 0; i < P; ++i) {
        pset[idx][0] = i;
        pset[idx][1] = j;
        ++idx;
      }
    }
  }
  else if constexpr (DIM == 3) {
    // 3D case: ix-fast, then iy, then iz
    // Loop order: k (outer), j (middle), i (inner) so i varies fastest
    int idx = 0;
    for (int k = 0; k < P; ++k) {
      for (int j = 0; j < P; ++j) {
        for (int i = 0; i < P; ++i) {
          pset[idx][0] = i;
          pset[idx][1] = j;
          pset[idx][2] = k;
          ++idx;
        }
      }
    }
  }

  return pset;
}

/**
 * @brief Evaluates the dir-th partial derivative of a tensor basis function.
 *
 * Product rule: d/dx_dir phi_p = phi_1d'(x_dir) * prod_{d!=dir} phi_1d(x_d),
 * with the 1D Legendre factors taken from pset[p].
 */
template <unsigned int DIM>
inline double
eval_tensor_basis_gradient (const std::array<double, DIM> &x, int p, int dir, const std::vector<multiindex<DIM>> &pset)
{
  double result = 1.0;
  for (unsigned int d = 0; d < DIM; ++d)
    result *= (d == static_cast<unsigned int> (dir)) ? phi_prime_1d (x[d], pset[p][d]) : phi_1d (x[d], pset[p][d]);

  return result;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

#pragma once

#ifdef T8_ENABLE_MRA

#include <array>
#include <vector>

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
    // 2D case: nested loops
    int idx = 0;
    for (int i = 0; i < P; ++i) {
      for (int j = 0; j < P; ++j) {
        pset[idx][0] = i;
        pset[idx][1] = j;
        ++idx;
      }
    }
  }
  else if constexpr (DIM == 3) {
    // 3D case: triple nested loops
    int idx = 0;
    for (int i = 0; i < P; ++i) {
      for (int j = 0; j < P; ++j) {
        for (int k = 0; k < P; ++k) {
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
 * @brief Evaluates a single tensor basis function at a point
 *
 * Computes phi_p(x) = product_{d=0}^{DIM-1} phi_1d(x[d], pset[p][d])
 * where phi_1d is the 1D Legendre basis function.
 *
 * @tparam DIM Spatial dimension
 * @param x Point in [0,1]^DIM where to evaluate
 * @param p Basis function index
 * @param pset Multiindex set mapping p to 1D indices
 * @param phi_1d Function pointer to 1D basis evaluation
 * @return double Value of the p-th basis function at x
 */
template <unsigned int DIM>
inline double
eval_tensor_basis (const std::array<double, DIM> &x, int p, const std::vector<multiindex<DIM>> &pset,
                   double (*phi_1d) (double, int))
{
  double result = 1.0;
  for (unsigned int d = 0; d < DIM; ++d)
    result *= phi_1d (x[d], pset[p][d]);

  return result;
}

/**
 * @brief Evaluates the gradient of a tensor basis function at a point
 *
 * Computes grad_phi_p(x)[dir] using the product rule:
 * d/dx_dir [phi_p(x)] = phi_1d'(x[dir], pset[p][dir]) * product_{d!=dir} phi_1d(x[d], pset[p][d])
 *
 * @tparam DIM Spatial dimension
 * @param x Point in [0,1]^DIM where to evaluate
 * @param p Basis function index
 * @param dir Derivative direction (0, 1, ..., DIM-1)
 * @param pset Multiindex set mapping p to 1D indices
 * @param phi_1d Function pointer to 1D basis evaluation
 * @param phi_prime_1d Function pointer to 1D basis derivative evaluation
 * @return double Value of the dir-th partial derivative of phi_p at x
 */
template <unsigned int DIM>
inline double
eval_tensor_basis_gradient (const std::array<double, DIM> &x, int p, int dir, const std::vector<multiindex<DIM>> &pset,
                            double (*phi_1d) (double, int), double (*phi_prime_1d) (double, int))
{
  double result = 1.0;
  for (unsigned int d = 0; d < DIM; ++d) {
    if (d == static_cast<unsigned int> (dir)) {
      // Use derivative in this direction
      result *= phi_prime_1d (x[d], pset[p][d]);
    }
    else {
      // Use function value in other directions
      result *= phi_1d (x[d], pset[p][d]);
    }
  }
  return result;
}

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <vector>
#include <stdexcept>

#include <t8_mra/num/vec.hpp>

namespace t8_mra
{

class mat {
  std::vector<double> data;
  size_t num_rows = 0u;
  size_t num_cols = 0u;

 public:
  mat () = default;
  mat (size_t _rows, size_t _cols): data (_rows * _cols, {}), num_rows (_rows), num_cols (_cols)
  {
  }

  mat (size_t _rows, size_t _cols, std::initializer_list<double> l): data (l), num_rows (_rows), num_cols (_cols)
  {
    if (l.size () != _rows * _cols)
      throw std::out_of_range ("number elements in t8_mra::util::mat does not fit to number columns "
                               "and number rows");
  }

  mat (const mat&) = default;
  mat&
  operator= (const mat&)
    = default;
  mat (mat&&) = default;
  mat&
  operator= (mat&&)
    = default;

  mat&
  operator= (const std::initializer_list<double>& l)
  {
    if (l.size () != num_rows * num_cols)
      throw std::out_of_range ("number elements in t8_mra::util::mat does not fit to number columns "
                               "and number rows");

    std::copy_n (l.begin (), l.size (), data.begin ());

    return *this;
  }

  double&
  operator() (size_t i, size_t j);
  double
  operator() (size_t i, size_t j) const;

  mat&
  operator= (double v);

  void
  resize (size_t _m, size_t _n);

  size_t
  rows () const noexcept;
  size_t
  cols () const noexcept;
};

inline double&
mat::operator() (size_t i, size_t j)
{
  if (i >= num_rows || j >= num_cols)
    throw std::out_of_range ("indices in t8_mra::util::mat::operator() is out of range");

  return data[num_rows * j + i];
}

inline double
mat::operator() (size_t i, size_t j) const
{
  if (i >= num_rows || j >= num_cols)
    throw std::out_of_range ("indices in t8_mra::util::mat::operator() is out of range");

  return data[num_rows * j + i];
}

inline void
mat::resize (size_t _rows, size_t _cols)
{
  data.clear ();
  num_rows = _rows;
  num_cols = _cols;
  data.resize (_rows * _cols);
}

inline size_t
mat::rows () const noexcept
{
  return num_rows;
}
inline size_t
mat::cols () const noexcept
{
  return num_cols;
}

/// Matrix is saved as A = (L - E_n) + U
/// below diagonal: L
/// Remaining matrix: U
inline void
lu_factors (mat& A, std::vector<size_t>& p)
{
  if (A.rows () != A.cols ())
    throw std::logic_error ("Matrix in t8_mra::util::lr_factor is not a square matrix");

  const auto n = A.rows ();
  p.resize (n);

  for (auto i = 0u; i < n; ++i)
    p[i] = i;

  for (auto j = 0u; j < n; j++) {
    auto Aj_max = 0.0;
    auto piv = j;

    for (auto k = j; k < n; k++) {
      const auto Ap = std::abs (A (k, j));

      if (Ap > Aj_max) {
        Aj_max = Ap;
        piv = k;
      }
    }

    if (piv != j) {
      std::swap (p[piv], p[j]);
      for (auto k = 0u; k < n; k++)
        std::swap (A (piv, k), A (j, k));
    }

    for (auto i = j + 1; i < n; i++) {
      A (i, j) /= A (j, j);

      for (auto k = j + 1; k < n; k++)
        A (i, k) -= A (i, j) * A (j, k);
    }
  }
}

inline void
lu_solve (const mat& A, const std::vector<size_t>& p, vec& x)
{
  if (A.rows () != A.cols ())
    throw std::logic_error ("Matrix in t8_mra::util::lr_solve is not a square matrix");
  if (A.rows () != p.size ())
    throw std::logic_error ("Permutation vector in t8_mra::util::lr_solve does not fit");
  if (A.rows () != x.size ())
    throw std::logic_error ("Solution vector in t8_mra::util::lr_solve does not fit");

  const auto n = A.rows ();

  const auto b = x;
  for (auto i = 0u; i < n; ++i) {
    x (i) = b (p[i]);
    for (auto k = 0u; k < i; ++k)
      x (i) -= A (i, k) * x (k);
  }

  for (int i = n - 1; i >= 0; --i) {
    for (auto k = static_cast<size_t> (i + 1); k < n; ++k)
      x (i) -= A (i, k) * x (k);
    x (i) /= A (i, i);
  }
}

}  // namespace t8_mra

#endif

#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

#include <cstddef>

namespace t8_mra
{

// ============================================================================
// Adding a new element shape
// ============================================================================
// Specialize shape_traits<Shape> below with the shape's compile-time mesh facts
// (dimension, child count, vertex count, VTK cell type, DOF count). The basis
// function space lives in basis<Shape, P> (num/basis/basis.hxx); the lmi bit layout
// in lmi_properties (data/levelmultiindex.hxx); projection and mask
// coefficients in the per-shape multiscale<> specialization (shapes/).

template <t8_eclass TShape>
concept is_cartesian = (TShape == T8_ECLASS_LINE || TShape == T8_ECLASS_QUAD || TShape == T8_ECLASS_HEX);

/// Binomial coefficient (compile-time), for simplex DOF counts.
constexpr size_t
binom (int n, int k) noexcept
{
  return (k > n)                  ? 0
         : (k == 0 || k == n)     ? 1
         : (k == 1 || k == n - 1) ? n
         : (2 * k < n)            ? binom (n - 1, k - 1) * n / k
                                  : binom (n - 1, k) * n / (n - k);
}

/// Per-shape compile-time facts. Primary template left undefined-ish (DIM 0)
/// so an unsupported shape fails loudly where the values are used.
template <t8_eclass TShape>
struct shape_traits
{
  static constexpr unsigned short DIM = 0;
  static constexpr unsigned short NUM_CHILDREN = 0;
  static constexpr int NUM_VERTICES = 0;
  static constexpr int VTK_CELL_TYPE = 0;

  static constexpr unsigned short
  dof (unsigned short)
  {
    return 0;
  }
};

template <>
struct shape_traits<T8_ECLASS_LINE>
{
  static constexpr unsigned short DIM = 1;
  static constexpr unsigned short NUM_CHILDREN = 2;
  static constexpr int NUM_VERTICES = 2;
  static constexpr int VTK_CELL_TYPE = 68;  // VTK_LAGRANGE_CURVE

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return P;
  }
};

template <>
struct shape_traits<T8_ECLASS_TRIANGLE>
{
  static constexpr unsigned short DIM = 2;
  static constexpr unsigned short NUM_CHILDREN = 4;
  static constexpr int NUM_VERTICES = 3;
  static constexpr int VTK_CELL_TYPE = 69;  // VTK_LAGRANGE_TRIANGLE

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return binom (DIM + P - 1, DIM);
  }
};

template <>
struct shape_traits<T8_ECLASS_QUAD>
{
  static constexpr unsigned short DIM = 2;
  static constexpr unsigned short NUM_CHILDREN = 4;
  static constexpr int NUM_VERTICES = 4;
  static constexpr int VTK_CELL_TYPE = 70;  // VTK_LAGRANGE_QUADRILATERAL

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return P * P;
  }
};

template <>
struct shape_traits<T8_ECLASS_HEX>
{
  static constexpr unsigned short DIM = 3;
  static constexpr unsigned short NUM_CHILDREN = 8;
  static constexpr int NUM_VERTICES = 8;
  static constexpr int VTK_CELL_TYPE = 72;  // VTK_LAGRANGE_HEXAHEDRON

  static constexpr unsigned short
  dof (unsigned short P)
  {
    return P * P * P;
  }
};

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

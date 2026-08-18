#pragma once

#ifdef T8_ENABLE_MRA

#include <algorithm>
#include <cstddef>

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

// ============================================================================
// Adding a new element shape
// ============================================================================
// Every per-shape specialization lives in one file per shape: shape_traits +
// mst policies in core/shape/, lmi layout in data/shape/, basis + cell_geometry
// in num/shape/, DG numerics in dg/. This header holds only the primary
// template and pulls the specializations in at the bottom.

template <t8_eclass TShape>
concept is_cartesian = (TShape == T8_ECLASS_LINE || TShape == T8_ECLASS_QUAD || TShape == T8_ECLASS_HEX);

/// Binomial coefficient (compile-time), for simplex DOF counts.
[[nodiscard]] constexpr size_t
binom (size_t n, size_t k) noexcept
{
  if (k > n)
    return 0;
  k = std::min (n - k, k);

  auto result = 1u;
  for (size_t i = 0; i < k; ++i)
    result = result * (n - i) / (i + 1);

  return result;
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

  [[nodiscard]] static constexpr unsigned short
  dof (unsigned short /*unused*/)
  {
    return 0;
  }
};

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/core/shape/cartesian.hxx"
#include "t8_mra/core/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA

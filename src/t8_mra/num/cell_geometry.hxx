#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

/// Reference-coordinate containment slack (boundary points, affine round-off).
inline constexpr double reference_cell_tol = 1e-9;

/// Cached affine geometry of one leaf. One specialization per shape (num/shape/);
/// trivially copyable.
template <t8_eclass Shape, int P>
struct cell_geometry;

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/num/shape/cartesian.hxx"
#include "t8_mra/num/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA

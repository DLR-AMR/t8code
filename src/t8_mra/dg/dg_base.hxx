#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

/// Per-shape DG numerics of a leaf: projection, evaluation and cell geometry.
/// t8code-free (operates on corner coordinates + volume). Specialized per shape
/// in dg/shape/.
template <t8_eclass TShape, int U, int P>
class dg;

}  // namespace t8_mra

// Per-shape specializations (defined after the primary template).
#include "t8_mra/dg/shape/cartesian.hxx"
#include "t8_mra/dg/shape/prism.hxx"
#include "t8_mra/dg/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA

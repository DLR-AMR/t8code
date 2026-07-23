#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

/// Per-shape DG numerics of a leaf: projection, evaluation and cell geometry.
/// t8code-free (operates on corner coordinates + volume). Specialized per shape
/// in dg/.
template <t8_eclass Shape, int U, int P>
class dg;

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

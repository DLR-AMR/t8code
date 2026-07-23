#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

/// Per-shape vertex-order handling under refinement. Specialized per shape in
/// core/shape/.
template <t8_eclass Shape>
struct ordering_policy;

/// Per-shape MST normalization and detail-norm scaling. Specialized per shape in
/// core/shape/.
template <t8_eclass Shape>
struct mst_scaling_policy;

}  // namespace t8_mra

#endif  // T8_ENABLE_MRA

#pragma once

#ifdef T8_ENABLE_MRA

#include <t8_eclass/t8_eclass.h>

namespace t8_mra
{

// ============================================================================
// Adding a new element shape
// ============================================================================
// Specialize vtk_shape<TShape> in io/shape/ and the writer in io/vtk.hxx
//
//   MAX_LAGRANGE_ORDER           highest Lagrange order
//   lagrange_nodes (order)       reference nodes in VTK's ordering
//   vertex_slot (corner, order)  VTK corner slot fed by t8code corner `corner`;
//   to_physical (ref, vertices)  geometry map over the VTK-ordered corners

/// Per-shape VTK node layout and geometry map. Specialized in io/shape/.
template <t8_eclass TShape>
struct vtk_shape;

}  // namespace t8_mra

/// Per-shape specializations (defined after the primary template).
#include "t8_mra/io/shape/cartesian.hxx"
#include "t8_mra/io/shape/prism.hxx"
#include "t8_mra/io/shape/triangle.hxx"

#endif  // T8_ENABLE_MRA

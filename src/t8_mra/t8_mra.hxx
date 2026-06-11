#pragma once

/** \file t8_mra.hxx
 * Main header of the t8code multiresolution analysis (MRA) module.
 *
 * Provides Harten-style adaptive grids based on multiscale transformations
 * of DG data:
 *
 *   t8_mra::multiscale<TShape, U, P>  with  TShape element shape,
 *                                           U solution components,
 *                                           P polynomial order
 *
 * Workflow: initialize_data() projects a function onto a uniform forest,
 * coarsen()/refine() adapt the grid driven by exchangeable criteria
 * (see criteria/), io/vtk.hxx writes the results.
 */

#ifdef T8_ENABLE_MRA

#include <t8_mra/shapes/triangle.hxx>
#include <t8_mra/shapes/cartesian.hxx>
#include <t8_mra/criteria/coarsening_criterion.hxx>
#include <t8_mra/criteria/refinement_criterion.hxx>
#include <t8_mra/io/vtk.hxx>

#endif  // T8_ENABLE_MRA

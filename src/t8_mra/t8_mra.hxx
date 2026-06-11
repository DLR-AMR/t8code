#pragma once

/** \file t8_mra.hxx
 * Main header of the t8code multiresolution analysis (MRA) module.
 *
 * The module maintains element-local DG data together with an adaptive
 * forest. A multiscale transformation decomposes the data into coarse
 * approximations and detail coefficients; adaptation criteria act on those
 * details: coarsening removes cells that carry no essential information,
 * refinement adds resolution where a criterion demands it. The criteria are
 * exchangeable (see criteria/), defaults implement Harten-style hard
 * thresholding and prediction.
 *
 * Central class: t8_mra::multiscale<TShape, U, P> with element shape TShape,
 * U solution components and polynomial order P.
 *
 * Basic usage:
 * \code
 * #include <t8_mra/t8_mra.hxx>
 *
 * t8_mra::multiscale<T8_ECLASS_TRIANGLE, U, P> mra (max_level, c_thresh, gamma, quadrature_rule, balanced, comm);
 *
 * // Project a function onto a uniform forest of level max_level
 * mra.initialize_data (cmesh, scheme, max_level, func);
 *
 * // Adapt: defaults are hard thresholding / Harten's prediction ...
 * mra.coarsen (min_level, max_level);
 * mra.refine (min_level, max_level);
 *
 * // ... or any type satisfying the coarsening/refinement_criterion concept
 * mra.refine (min_level, max_level, my_criterion {});
 *
 * t8_mra::write_forest_lagrange_vtk (mra, "solution", P - 1);
 *
 * mra.cleanup ();
 * \endcode
 */

#ifdef T8_ENABLE_MRA

#include <t8_mra/shapes/triangle.hxx>
#include <t8_mra/shapes/cartesian.hxx>
#include <t8_mra/criteria/coarsening_criterion.hxx>
#include <t8_mra/criteria/refinement_criterion.hxx>
#include <t8_mra/io/vtk.hxx>

#endif  // T8_ENABLE_MRA

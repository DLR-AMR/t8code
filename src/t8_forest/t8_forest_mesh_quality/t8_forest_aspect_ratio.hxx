/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2026 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file t8_forest_aspect_ratio.hxx
 * MEasure the aspect ratio of the elements of a forest.
 *
 * The aspect ratio (AR) measures how much an element differs from the ideal (regular) element shape.
 *
 *  AR = L_ax * |boundary(E)| / (d * |E|) / c_shape,
 *
 * where L_max is the longest edge of the element E, |boundary(E)| the perimeter (2D) or the surface area (3D), |E| it area (2D) or volume(3D),
 * and d the dimension.
 * The constant c_shape is chosen so that the ideal element has AR = 1.
 * d * |E| / |boundary(E)| is the inradius of triangles, squares, tetrahedra and cubes, so AR is the longest edge divided by the inradius,
 * normalized to 1.
 *
 * Ideal elements (AR = 1):
 *  - triangle: equilateral triangle,
 *  - quad:     square,
 *  - tet:      regular tetrahedron,
 *  - hex:      cube.
 *  - prism:    right prism over an equilateral triangle whose edges all have the same length,
 *  - pyramid:  square pyramid whose edges all have the same length.
 *
 * The more distorted (stretched, sheared or flattened) an element is, the larger its AR.
 * A degenerate element (zero area/volume) has AR = infinity.
 * AR >= 1 holds for triangles, quads, tets, hexes and prisms. Pyramids are an exception: a pyramid that is taller than the ideal one can have AR < 1.
 * For straight-sided triangles, flat quads and tetrahedra this is the same definition as the "aspect ratio" of the Verdict library,
 * which ParaView's "Mesh Quality" filter uses. Values can therefore be cross-checked there.
 *
 * Some other tools report a quality only with [0, 1] instead, where 1 is ideal and 0 is degenerate. In this file that is 1/AR, which is what the histogram
 * of \ref t8_forest_compute_aspect_ratio_stats and \ref t8_forest_print_aspect_ratio uses.
 *
 * The aspect ratio is computed from the corner coordinates of the element.
 * Curved geometries are therefore treated as if the element had straight edges between the corner points.
 * Faces that are not flat (warped quad faces) are approximated by the average of their two possible splits into two triangles.
 *
 * In the hxx, we always declare the functions and implement them in the cxx.
 * In some files we also provide a c_interface or .c files. We do not provide them here as we only want
 * to ensure that this code is callable from C++ and not also for C.
 */

#pragma once

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <array>

/** Number of bins of the aspect ratio histogram in \ref t8_forest_aspect_ratio_stats
 *  The bins divide the quality 1/AR in [0, 1] into equally sized intervals,
 *  as in gmsh's quality output: 0.0 <= 1/AR <0.1, 0.1 <= 1/AR <0.2, and so on.
 */
inline constexpr int T8_FOREST_ASPECT_RATIO_HISTOGRAM_BINS = 10;

/** Statistics of the aspect ratios of all leaf elements of a forest.
 *  See \ref t8_forest_compute_aspect_ratio_stats.
 */
struct t8_forest_aspect_ratio_stats
{
  t8_gloidx_t num_elements; /**< Global number of leaf elements in the forest. */
  double min;               /**< Smallest aspect ratio, i.e. the AR of best element (0 if forest is empty) */
  double max;               /**< Largest aspect ratio, i.e. the AR of worst element (0 if forest is empty) */
  double average;           /**< Arithmetic mean of the aspect ratios (0 if the forest is empty) */
  /** histogram[i] is the global number of elements whose quality 1/AR lies in
   *  [i / N, (i + 1) / N)], with N = \ref T8_FOREST_ASPECT_RATIO_NUM_BINS.
   *  The last bin also contains the elements with 1/AR >= 1 (ideal elements, and pyramids with AR slightly below 1).
   *  The sum of all bins is \a num_elements.
   */
  std::array<t8 : gloidx_t, T8_FOREST_ASPECT_RATIO_NUM_BINS> histogram;
};

/** Compute the aspect ratio of a single element.
 *  See the file description for a definition. Vertices and lines cannot be distorted, so their aspect ratio is 1.
 *  \param [in] forest  The forest. Must be committed.
 *  \param [in] ltreeid The forest local id of the tree in which the element is.
 *  \param [in] element The element.
 *  \return             The aspect ratio of \a element. It is 1 for an ideal element an >= 1 for all shapes except pyramids (see the file description).
 *                      Larger values mean a more distorted element. It is infinity for a degenerate element.
 */
double
t8_forest_element_aspect_ratio (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/** Compute statistics (minimum, maximum, average and a histogram) of the aspect ratios of all leaf elements of a forest.
 *  This function is collective: it must be called on all processes of the forest's communicator, and every process receives the global result.
 *  \param [in] forest  The forest. Must be committed.
 *  \return             The global aspect ratio statistics of \a forest.
 */
t8_forest_aspect_ratio_stats
t8_forest_compute_aspect_ratio_stats (t8_forest_t forest);

/** Compute the average aspect ratio of all leaf elements of a forest.
 *  This function is collective
 *  \param [in] forest  The forest. Must be committed.
 *  \return             The arithmetic mean of the aspect ratios of all leaf elements of \a forest, or 0 if the forest has no elements.
 */
double
t8_forest_average_aspect_ratio (t8_forest_t forest);

/** Compute the largest aspect ratio of all leaf elements of a forest, i.e. the aspect ratio of the worst element of the mesh.
 *  This function is collective.
 *  \param [in] forest  The forest. Must be committed.
 *  \return             The maximum aspect ratio of all leaf elements of \a forest, or 0 if the forest has no elements.
 */
double
t8_forest_max_aspect_ratio (t8_forest_t forest);

/** Print an overview of the aspect ratios of all leaf elements of a forest:
 *  the best, worst and average aspect ratio, and a gmsh-style histogram that shows
 *  how many elements (and which percentage) lie in each quality range 1/AR.
 *  This function is collective. The output is only printed on rank 0.
 * \param [in] forest   The forest. Must be committed.
 */
void
t8_forest_print_aspect_ratio (t8_forest_t forest);

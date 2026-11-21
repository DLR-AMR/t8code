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

/** \file t8_forest_aspect_ratio.cxx
 *  Implementation of the aspect ratio functions declared in \ref t8_forest_aspect_ratio.hxx.
 */

#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_eclass/t8_eclass.h>
#include "t8_forest_aspect_ratio.hxx"
#include <t8_schemes/t8_scheme.hxx>
#include <t8_forest/t8_forest_geometrical.h> /* For t8_forest_element_coordinate. */
#include <t8_vector/t8_vec.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>

/** Compute L_max * |boundary(E)| / (d * |E|) for the ideal element E of a given shape,
 *  with L_max the longest edge of the ideal element, |boundary(E)| the measure of the element boundary,
 *  d the spatial dimension, and |E| the measure of the element volume.
 *  We divide by this value so that the ideal element has an aspect ratio of 1.
 *  \param [in] shape The shape of the element. Mus be 2D or 3D.
 *  \return           The value for the ideal element of \a shape.
 */
static double
t8_forest_aspect_ratio_ideal_value (const t8_element_shape_t shape)
{
  switch (shape) {
  case T8_ECLASS_TRIANGLE:
    /* Equilateral triangle with edge length 1: perimeter 3, area sqrt(3)/4. */
    return 2.0 * std::sqrt (3.0);
  case T8_ECLASS_QUAD:
    /* Unit square: perimeter 4, area 1. */
    return 2.0;
  case T8_ECLASS_TET:
    /* Regular tetrahedron with edge length 1: surface area sqrt (3), volume 1/(6*sqrt(2)). */
    return 2.0 * std::sqrt (6.0);
  case T8_ECLASS_HEX:
    /* Unit cube: surface area 6, volume 1. */
    return 2.0;
  case T8_ECLASS_PRISM:
    /* Prism with all edges of length 1: surface 3 + sqrt(3)/2, volume sqrt(3)/4. */
    return 2.0 / 3.0 + 4.0 / std::sqrt (3.0);
  case T8_ECLASS_PYRAMID:
    /* Square pyramid with all edges of length 1: surface 1 + sqrt(3), volume 1/(3*sqrt(2)). */
    return std::sqrt (2.0) * (1.0 + std::sqrt (3.0));
  default:
    SC_ABORT_NOT_REACHED ();
  }
  return 0.0;
}

/** Compute the area of a triangle.
 *  \param [in] a, b, c The corners of the triangle.
 *  \return             The area of the triangle.
 */
static double
t8_forest_aspect_ratio_triangle_area (const t8_3D_vec &a, const t8_3D_vec &b, const t8_3D_vec &c)
{
  t8_3D_vec ab, ac, cross;
  t8_diff (b, a, ab);
  t8_diff (c, a, ac);
  t8_cross_3D (ab, ac, cross);
  return 0.5 * t8_norm (cross);
}

/** Compute the volume of a tetrahedron.
 *  \param [in] a, b, c, d The corners of the tetrahedron.
 *  \return                 The volume of the tetrahedron.
 */
static double
t8_forest_aspect_ratio_tet_volume (const t8_3D_vec &a, const t8_3D_vec &b, const t8_3D_vec &c, const t8_3D_vec &d)
{
  t8_3D_vec ab, ac, ad, cross;
  t8_diff (b, a, ab);
  t8_diff (c, a, ac);
  t8_diff (d, a, ad);
  t8_cross_3D (ab, ac, cross);
  return std::fabs (t8_dot (cross, ad)) / 6.0;
}

/** Compute the histogram bin of an aspect ratio, see \ref
 *  t8_forest_aspect_ratio_stat::histogram.
 *  \param [in] aspect_ratio The aspect ratio.
 *  \return                  The bin in [0, T8_FOREST_ASPECT_RATIO_NUM_BINS - 1] that \a aspect_ratio belongs to.
 */
static int
t8_forest_aspect_ratio_bin (const double aspect_ratio)
{
  T8_ASSERT (!std::isnan (aspect_ratio));
  constexpr int num_bins = T8_FOREST_ASPECT_RATIO_NUM_BINS;
  /* The quality 1/AR lies in [0, 1] (0 for degenerate elements with AR = infinity).
   * It can be slightly above 1 because of rounding errors, or for pyramids (see t8_forest_aspect_ratio.hxx).
   * We clamp before casting to int, so that these elements end up in the last bin.
   */
  const double scaled_quality = std::clamp (num_bins / aspect_ratio, 0.0, num_bins - 1.0);
  return static_cast<int> (scaled_quality);
}

double
t8_forest_element_aspect_ratio (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element)
{
  T8_ASSERT (t8_forest_is_committed (forest));
  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  T8_ASSERT (scheme->element_is_valid (tree_class, element));

  /* The shape of an element can differ from the class of its tree (a pyramid tree also contains tets), so we ask the scheme. */
  const t8_element_shape_t shape = scheme->element_get_shape (tree_class, element);
  const int dim = t8_eclass_to_dimension[shape];
  if (dim < 2) {
    /* Vertices and lines cannot be distorted. */
    return 1.0;
  }

  /* Compute the corner coordinates and the center (mean of the corners) of the element. */
  const int num_corners = scheme->element_num_corners (tree_class, element);
  T8_ASSERT (0 < num_corners && num_corners <= T8_ECLASS_MAX_CORNERS);
  std::array<t8_3D_vec, T8_ECLASS_MAX_CORNERS> corners;
  t8_3D_vec center = { 0.0, 0.0, 0.0 };
  for (int icorner = 0; icorner < num_corners; ++icorner) {
    t8_forest_element_coordinate (forest, ltreeid, element, icorner, corners[icorner].data ());
    t8_axpy (corners[icorner], center, 1.0);
  }
  t8_ax (center, 1.0 / num_corners);

  /* We loop over the faces of the element (the edges in 2D) to compute
   *  - the length of the longest edge. Every edge of the element lies on a face.
   *    In 3D each edge lies on two faces and is visited twice, which does not change the maximum.
   *  - the measure of the boundary (perimeter in 2D, surface area in 3D).
   *  - the measure of the element (area in 2D, volume in 3D). We decompose the element into
   *    one triangle (2D) or tetrahedron (3D) per face (quad faces: two tetrahedra) with the center as apex.
   *    This is exact for convex elements with flat faces.
   *    We ask the scheme for the corners of each face, since the face numbering can depend on the element type. */
  double longest_edge = 0.0;
  double boundary = 0.0;
  double measure = 0.0;
  const int num_faces = scheme->element_get_num_faces (tree_class, element);
  for (int iface = 0; iface < num_faces; ++iface) {
    const t8_element_shape_t face_shape = scheme->element_get_face_shape (tree_class, element, iface);
    const int num_face_corners = t8_eclass_num_vertices[face_shape];
    T8_ASSERT (num_face_corners <= T8_ECLASS_MAX_CORNERS_2D);
    std::array<int, T8_ECLASS_MAX_CORNERS_2D> face_corners;
    for (int iface_corner = 0; iface_corner < num_face_corners; ++iface_corner) {
      face_corners[iface_corner] = scheme->element_get_face_corner (tree_class, element, iface, iface_corner);
      T8_ASSERT (0 <= face_corners[iface_corner] && face_corners[iface_corner] < num_corners);
    }

    switch (face_shape) {
    case T8_ECLASS_LINE: {
      /* The faces of a 2D element are its edges. */
      const t8_3D_vec &a = corners[face_corners[0]];
      const t8_3D_vec &b = corners[face_corners[1]];
      const double length = t8_dist (a, b);
      longest_edge = std::max (longest_edge, length);
      boundary += length;
      measure += t8_forest_aspect_ratio_triangle_area (center, a, b);
      break;
    }
    case T8_ECLASS_TRIANGLE: {
      const t8_3D_vec &a = corners[face_corners[0]];
      const t8_3D_vec &b = corners[face_corners[1]];
      const t8_3D_vec &c = corners[face_corners[2]];
      longest_edge = std::max ({ longest_edgem t8_dist (a, b), t8_dist (b, c), t8_dist (c, a) });
      boundary += t8_forest_aspect_ratio_triangle_area (a, b, c);
      measure += t8_forest_aspect_ratio_tet_volume (center, a, b, c);
      break;
    }
    case T8_ECLASS_QUAD: {
      const t8_3D_vec &a = corners[face_corners[0]];
      const t8_3D_vec &b = corners[face_corners[1]];
      const t8_3D_vec &c = corners[face_corners[2]];
      const t8_3D_vec &d = corners[face_corners[3]];
      longest_edge = std::max ({ longest_edge, t8_dist (a, b), t8_dist (b, d), t8_dist (d, c), t8_dist (c, a) });
      const double area_split_ad
        = t8_forest_aspect_ratio_triangle_area (a, b, d) + t8_forest_aspect_ratio_triangle_area (a, d, c);
      const double area_split_bc
        = t8_forest_aspect_ratio_triangle_area (a, b, c) + t8_forest_aspect_ratio_triangle_area (b, d, c);
      boundary += 0.5 * (area_split_ad + area_split_bc);
      const double volume_split_ad
        = t8_forest_aspect_ratio_tet_volume (center, a, b, d) + t8_forest_aspect_ratio_tet_volume (center, a, d, c);
      const double volume_split_bc
        = t8_forest_aspect_ratio_tet_volume (center, a, b, c) + t8_forest_aspect_ratio_tet_volume (center, b, d, c);
      measure += 0.5 * (volume_split_ad + volume_split_bc);
      break;
    }
    default:
      SC_ABORT_NOT_REACHED ();
    }
  }

  if (!(measure > 0.0)) {
    /* The element is degenerate (zero area or volume). */
    return std::numeric_limits<double>::infinity ();
  }
  return longest_edge * boundary / (dim * measure * t8_forest_aspect_ratio_ideal_value (shape));
}

t8_forest_aspect_ratio_stats
t8_forest_compute_aspect_ratio_stats (t8_forest_t forest)
{
  T8_ASSERT (t8_forest_is_committed (forest));

  /* Compute the statistics of the process local elements. */
  double local_sum = 0.0;
  double local_min = std::numeric_limits<double>::infinity ();
  double local_max = -std::numeric_limits<double>::infinity ();
  std::array<t8_gloidx_t, T8_FOREST_ASPECT_RATIO_NUM_BINS> local_histogram = {};

  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);
  for (t8_locidx_t itree = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement) {
      const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);
      const double aspect_ratio = t8_forest_element_aspect_ratio (forest, itree, element);
      local_sum += aspect_ratio;
      local_min = std::min (local_min, aspect_ratio);
      local_max = std::max (local_max, aspect_ratio);
      ++local_histogram[t8_forest_aspect_ratio_bin (aspect_ratio)];
    }
  }

  /* Combine the results of all processes. */
  const sc_MPI_Comm comm = t8_forest_get_mpicomm (forest);
  t8_forest_aspect_ratio_stats stats;
  double global_sum;
  int mpiret = sc_MPI_Allreduce (&local_sum, &global_sum, 1, sc_MPI_DOUBLE, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (&local_min, &stats.min, 1, sc_MPI_DOUBLE, sc_MPI_MIN, comm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (&local_max, &stats.max, 1, sc_MPI_DOUBLE, sc_MPI_MAX, comm);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Allreduce (local_histogram.data (), stats.histogram.data (), T8_FOREST_ASPECT_RATIO_NUM_BINS,
                             sc_MPI_UINT64_T, sc_MPI_SUM, comm);
  SC_CHECK_MPI (mpiret);

  stats.num_elements = 0;
  for (const t8_gloidx_t count : stats.histogram) {
    stats.num_elements += count;
  }
  T8_ASSERT (stats.num_elements == t8_forest_get_global_num_leaf_elements (forest));

  if (stats.num_elements > 0) {
    stats.average = global_sum / stats.num_elements;
  }
  else {
    stats.min = 0.0;
    stats.max = 0.0;
    stats.average = 0.0;
  }
  return stats;
}

double
t8_forest_average_aspect_ratio (t8_forest_t forest)
{
  return t8_forest_compute_aspect_ratio_stats (forest).average;
}

double
t8_forest_max_aspect_ratio (t8_forest_t forest)
{
  return t8_forest_compute_aspect_ratio_stats (forest).max;
}

void
t8_forest_print_aspect_ratio (t8_forest_t forest)
{
  const t8_forest_aspect_ratio_stats stats = t8_forest_compute_aspect_ratio_stats (forest);
  constexpr int num_bins = T8_FOREST_ASPECT_RATIO_NUM_BINS;

  t8_global_productionf ("Aspect ratio (AR) of the %" T8_GLOIDX_FORMAT
                         " leaf elements of the forest (AR = 1 is ideal, larger is worse):\n",
                         stats.num_elements);
  if (stats.num_elements == 0) {
    return;
  }
  t8_global_productionf (" best AR:    %.4f\n", stats.min);
  t8_global_productionf (" worst AR:   %.4f\n", stats.max);
  t8_global_productionf (" average AR: %.4f\n", stats.average);
  t8_global_productionf (" Distribution of the element quality 1/AR (1 = ideal, 0 = degenerate):\n");
  for (int ibin = 0; ibin < num_bins; ++ibin) {
    const double lower = static_cast<double> (ibin) / num_bins;
    const double upper = static_cast<double> (ibin + 1) / num_bins;
    const char *upper_relation = (ibin == num_bins - 1) ? "<=" : "<";

    char ar_range[64];
    if (ibin == 0) {
      snprintf (ar_range, sizeof (ar_range), "AR > %.2f", 1.0 / upper);
    }
    else if (ibin == num_bins - 1) {
      snprintf (ar_range, sizeof (ar_range), "AR <= %.2f", 1.0 / lower);
    }
    else {
      snprintf (ar_range, sizeof (ar_range), "%.2f < AR <= %.2f", 1.0 / upper, 1.0 / lower);
    }

    const double percentage = 100.0 * stats.histogram[ibin] / stats.num_elements;
    t8_global_productionf ("   %.2f <= 1/AR %s %.2f (%-18s): %10" T8_GLOIDX_FORMAT " elements (%6.2f %%)\n", lower,
                           upper_relation, upper, ar_range, stats.histogram[ibin], percentage);
  }
}

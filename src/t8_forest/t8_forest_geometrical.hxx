/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2025 the developers

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

/** \file t8_forest_geometrical.hxx
 * We define the geometrical queries for a forest of trees in this file.
 * Also includes \ref t8_forest_geometrical.h since not all forest functions are CPP.
 */

#ifndef T8_FOREST_GEOMETRICAL_HXX
#define T8_FOREST_GEOMETRICAL_HXX

#include <t8_forest/t8_forest_geometrical.h>
#include <t8_types/t8_vec.hxx>
#include <t8_vector_helper/t8_zip.hxx>
#include <t8_schemes/t8_scheme.hxx>

/** Compute the coordinates of a point inside an element inside a tree.
 *  The point is given in reference coordinates inside the element and gets
 *  converted to reference coordinates inside the tree. After that, the point
 *  is converted to global coordinates inside the domain. If needed, the element
 *  is stretched by the given stretch factors (the resulting mesh is then
 *  no longer non-overlapping).
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        A list of reference coordinates of points inside the element.
 * \param [out]     coords_out        An allocated list of x, y and z coordinates which will be filled.
 * \param [in]      stretch_factors   If provided, elements are stretched according to the stretch factors
 *                                    of the tree.
 * \tparam TRefCoords                 A container of t8_point with at least the same dimension as \a element.
 * \tparam TOutCoords                 A container of t8_point with dimension 3.
 */
template <T8PointContainerType TRefCoords, T8PointContainerType<3> TOutCoords>
void
t8_forest_element_from_ref_coords_ext (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                       const TRefCoords &ref_coords, TOutCoords &coords_out,
                                       const double *stretch_factors)
{
  T8_ASSERTF (std::ranges::size (ref_coords) <= std::ranges::size (coords_out),
              "coords_out is smaller than ref_coords.");

  const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, ltreeid);
  const int tree_dim = t8_eclass_to_dimension[tree_class];

  T8_ASSERTF (std::tuple_size_v<TRefCoords::value_type> >= tree_dim == 0 ? 1 : tree_dim,
              "Input reference coordinates do not have the right dimension.");

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  const t8_cmesh_t cmesh = t8_forest_get_cmesh (forest);
  const t8_gloidx_t gtreeid = t8_forest_global_tree_id (forest, ltreeid);

  /* Create an intermediate container with the same form as ref_coords */
  TRefCoords tree_ref_coords (ref_coords.size (), TRefCoords::value_type (ref_coords[0].size ()));

  if (stretch_factors != NULL) {
#if T8_ENABLE_DEBUG
    const t8_geometry_type_t geom_type = t8_geometry_get_type (cmesh, gtreeid);
    T8_ASSERT (geom_type == T8_GEOMETRY_TYPE_LINEAR || geom_type == T8_GEOMETRY_TYPE_LINEAR_AXIS_ALIGNED);
#endif /* T8_ENABLE_DEBUG */
    TRefCoords stretched_ref_coords (ref_coords.size (), TRefCoords::value_type (ref_coords[0].size ()));
    for (auto &[coord, stretched_coord] : std::ranges::views::zip (ref_coords, stretched_ref_coords)) {
      int8_t dim = 0;
      for (auto &[val, stretched_val] : std::ranges::views::zip (coord, stretched_coord)) {
        stretched_val = 0.5 + ((val - 0.5) * stretch_factors[dim]);
        ++dim;
      }
    }
    scheme->element_get_reference_coords (tree_class, element, stretched_ref_coords, tree_ref_coords);
  }
  else {
    scheme->element_get_reference_coords (tree_class, element, ref_coords, tree_ref_coords);
  }

  /* Conversion is only here to split the PR in two. The geometries will also be altered and then there will be no conversion */
  std::vector<double> converted_ref_coords;
  converted_ref_coords.reserve (tree_dim == 0 ? 1 : tree_dim);
  size_t count = 0;
  for (const auto &coord : tree_ref_coords) {
    for (int dim = 0; dim < tree_dim == 0 ? 1 : tree_dim; ++dim)
      converted_ref_coords.push_back (coord[dim]);
  }
  T8_ASSERT (converted_ref_coords.size () == std::ranges::size (ref_coords) * tree_dim == 0 ? 1 : tree_dim);
  /* End of conversion */

  t8_geometry_evaluate (cmesh, gtreeid, converted_ref_coords.data (), std::ranges::size (ref_coords), coords_out);

  T8_FREE (tree_ref_coords);
}

/** Compute the coordinates of a point inside an element inside a tree.
 *  The point is given in reference coordinates inside the element and gets
 *  converted to reference coordinates inside the tree. After that, the point
 *  is converted to global coordinates inside the domain.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        A list of reference coordinates of points inside the element.
 * \param [out]     coords_out        An allocated list of x, y and z coordinates which will be filled.
 * \tparam TRefCoords                 A container of t8_point with at least the same dimension as \a element.
 * \tparam TOutCoords                 A container of t8_point with dimension 3.
 */
template <T8PointContainerType TRefCoords, T8PointContainerType<3> TOutCoords>
void
t8_forest_element_from_ref_coords (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const TRefCoords &ref_coords, TOutCoords &coords_out)
{
  t8_forest_element_from_ref_coords_ext (forest, ltreeid, element, ref_coords, coords_out, NULL);
}

#endif /* T8_FOREST_GEOMETRICAL_HXX */

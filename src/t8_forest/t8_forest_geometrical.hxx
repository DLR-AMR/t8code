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

/** \file t8_forest_geometrical.hxx
 * C++ interface for the geometrical queries of a forest. Also includes \ref t8_forest_geometrical.h.
 */

#ifndef T8_FOREST_GEOMETRICAL_HXX
#define T8_FOREST_GEOMETRICAL_HXX

#include <t8_forest/t8_forest_geometrical.h>
#include <t8_types/t8_vec.hxx>
#include <span>

/** Compute the coordinates of a given vertex of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      corner_number     The corner number, in Z-order, of the vertex which should be computed.
 * \return                            The x, y and z coordinates of the vertex inside the domain.
 */
t8_3D_vec
t8_forest_element_coordinate (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                              const int corner_number);

/** Compute the coordinates of a point inside an element inside a tree.
 * Arrays with typical reference coordinates can be found in \ref t8_element.h.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        The reference coordinates of the point inside the element.
 * \return                            The coordinates of the point inside the domain.
 */
t8_3D_vec
t8_forest_element_from_ref_coords (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                   const t8_3D_vec &ref_coords);

/** Compute the coordinates of points inside an element inside a tree.
 * Arrays with typical reference coordinates can be found in \ref t8_element.h.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        The reference coordinates of the points inside the element.
 * \param [out]     coords_out        The coordinates of the points inside the domain.
 *                                    Must have the same size as \a ref_coords.
 */
void
t8_forest_element_from_ref_coords (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                   std::span<const t8_3D_vec> ref_coords, std::span<t8_3D_vec> coords_out);

/** Compute the coordinates of points inside an element inside a tree.
 * If needed, the element is stretched by the given stretch factors (the resulting mesh is then
 * no longer non-overlapping).
 * Arrays with typical reference coordinates can be found in \ref t8_element.h.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        The reference coordinates of the points inside the element.
 * \param [out]     coords_out        The coordinates of the points inside the domain.
 *                                    Must have the same size as \a ref_coords.
 * \param [in]      stretch_factors   If not nullptr, an array of at least d doubles (d = dimension of the tree).
 *                                    The element is stretched by these factors around the point with all
 *                                    reference coordinates equal to 0.5. For lines, quads and hexes, this is
 *                                    the center of the element.
 *                                    Only supported for linear geometries.
 */
void
t8_forest_element_from_ref_coords_ext (const t8_forest_t forest, const t8_locidx_t ltreeid, const t8_element_t *element,
                                       std::span<const t8_3D_vec> ref_coords, std::span<t8_3D_vec> coords_out,
                                       const double *stretch_factors);

#endif /* !T8_FOREST_GEOMETRICAL_HXX */

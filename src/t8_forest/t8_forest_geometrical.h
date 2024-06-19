/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

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

/** \file t8_forest_geometrical.h
 * We define the geometrical queries for a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest */

#ifndef T8_FOREST_GEOMETRICAL_H
#define T8_FOREST_GEOMETRICAL_H

#include <sc_statistics.h>
#include <t8_cmesh.h>
T8_EXTERN_C_BEGIN ();

/** Return the dimension of a forest.
 * \param [in]  forest    A forest.
 * \return                The dimension.
 * \a forest must be committed before calling this function.
 * Note: The dimension is inferred from the associated \b cmesh.
 */
int
t8_forest_get_dimension (const t8_forest_t forest);

/** Compute the coordinates of a given vertex of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      corner_number The corner number, in Z-order, of the vertex which should be computed.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the vertex.
 */
void
t8_forest_element_coordinate (t8_forest_t forest, t8_locidx_t ltree_id, const t8_element_t *element, int corner_number,
                              double *coordinates);

/** Compute the coordinates of a point inside an element inside a tree.
 *  The point is given in reference coordinates inside the element and gets
 *  converted to reference coordinates inside the tree. After that, the point
 *  is converted to global coordinates inside the domain. If needed, the element
 *  is stretched by the given stretch factors (the resulting mesh is then 
 *  no longer non-overlapping).
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        The reference coordinates of the point inside the element.
 * \param [in]      num_coords        The number of coordinate sets in ref_coord (dimension x double).
 * \param [out]     coords_out        On input an allocated array to store 3 doubles, on output
 *                                    the x, y and z coordinates of the point inside the domain.
 * \param [in]      stretch_factors   If provided, elements are stretched according to the stretch factors
 *                                    of the tree.
 */

void
t8_forest_element_from_ref_coords_ext (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                       const double *ref_coords, const size_t num_coords, double *coords_out,
                                       const double *stretch_factors);

/** Compute the coordinates of a point inside an element inside a tree.
 *  The point is given in reference coordinates inside the element and gets
 *  converted to reference coordinates inside the tree. After that, the point
 *  is converted to global coordinates inside the domain.
 * \param [in]      forest            The forest.
 * \param [in]      ltreeid           The forest local id of the tree in which the element is.
 * \param [in]      element           The element.
 * \param [in]      ref_coords        The reference coordinates of the point inside the element.
 * \param [in]      num_coords        The number of coordinate sets in ref_coord (dimension x double).
 * \param [out]     coords_out        On input an allocated array to store 3 doubles, on output
 *                                    the x, y and z coordinates of the point inside the domain.
 */
void
t8_forest_element_from_ref_coords (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element,
                                   const double *ref_coords, const size_t num_coords, double *coords_out);

/** Compute the coordinates of the centroid of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * The centroid can be seen as the midpoint of an element and thus can for example be used
 * to compute level-set values or the distance between two elements.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the centroid.
 */
void
t8_forest_element_centroid (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, double *coordinates);

/** Compute the diameter of an element if a geometry for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       For lines the value is exact while for other element types it is only
 *                             an approximation.
 */
double
t8_forest_element_diam (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/** Compute the volume of an element if a geometry for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       This function assumes d-linear interpolation for the
 *                             tree vertex coordinates.
 *                             \a forest must be committed when calling this function.
 */
double
t8_forest_element_volume (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element);

/** Compute the area of an element's face if a geometry for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \return                     The area of \a face.
 * \a forest must be committed when calling this function.
 */
double
t8_forest_element_face_area (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face);

/** Compute the vertex coordinates of the centroid of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the centroid of \a face.
 * \a forest must be committed when calling this function.
 */
void
t8_forest_element_face_centroid (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                                 double centroid[3]);

/** Compute the normal vector of an element's face if a geometry for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the normal vector of \a element at \a face.
 * \a forest must be committed when calling this function.
 */
void
t8_forest_element_face_normal (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, int face,
                               double normal[3]);

/** Compute whether a given face of a given leaf in a forest lies at the domain boundary.
 * This includes inner boundaries if the forest has deleted elements.
 * \param [in]      forest    The forest.
 * \param [in]      local_tree A local tree of \a forest.
 * \param [in]      leaf      A leaf element of \a local_tree in \a forest.
 * \param [in]      face      The face number of a face of \a leaf to check.
 * \return True (non-zero) if face \a face of \a leaf lies at the domain boundary.
*/
int
t8_forest_leaf_is_boundary (const t8_forest_t forest, t8_locidx_t local_tree, const t8_element_t *leaf, int face);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GEOMETRICAL_H */

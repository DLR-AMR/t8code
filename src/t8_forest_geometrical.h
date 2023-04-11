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

/*TODO: Check if all headers are needed */
#include <sc_statistics.h>
#include <t8_cmesh.h>
#include <t8_element.h>
#include <t8_vtk.h>
#include <t8_data/t8_containers.h>

/** This type controls, which neighbors count as ghost elements.
 * Currently, we support face-neighbors. Vertex and edge neighbors
 * will eventually be added. */
typedef enum
{
  T8_GHOST_NONE = 0,  /**< Do not create ghost layer. */
  T8_GHOST_FACES,     /**< Consider all face (codimension 1) neighbors. */
  T8_GHOST_EDGES,     /**< Consider all edge (codimension 2) and face neighbors. */
  T8_GHOST_VERTICES   /**< Consider all vertex (codimension 3) and edge and face neighbors. */
} t8_ghost_type_t;

T8_EXTERN_C_BEGIN ();

/** Compute the coordinates of a given vertex of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      corner_number The corner number, in Z-order, of the vertex which should be computed.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the vertex.
 */
void                t8_forest_element_coordinate (t8_forest_t forest,
                                                  t8_locidx_t ltree_id,
                                                  const t8_element_t *element,
                                                  int corner_number,
                                                  double *coordinates);

/** Compute the coordinates of the centroid of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * The centroid is the sum of all corner vertices divided by the number of corners.
 * The centroid can be seen as the midpoint of an element and thus can for example be used
 * to compute level-set values or the distance between two elements.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [out]     coordinates On input an allocated array to store 3 doubles, on output
 *                             the x, y and z coordinates of the centroid.
 */
void                t8_forest_element_centroid (t8_forest_t forest,
                                                t8_locidx_t ltreeid,
                                                const t8_element_t *element,
                                                double *coordinates);

/** Compute the diameter of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       For lines the value is exact while for other element types it is only
 *                             an approximation.
 */
double              t8_forest_element_diam (t8_forest_t forest,
                                            t8_locidx_t ltreeid,
                                            const t8_element_t *element);

/** Compute the volume of an element if a geometry
 * for this tree is registered in the forest's cmesh.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \return                     The diameter of the element.
 * \note                       This function assumes d-linear interpolation for the
 *                             tree vertex coordinates.
 *                             \a forest must be committed when calling this function.
 */
double              t8_forest_element_volume (t8_forest_t forest,
                                              t8_locidx_t ltreeid,
                                              const t8_element_t *element);

/** Compute the area of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * This is only an approximation.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \return                     The area of \a face.
 * \a forest must be committed when calling this function.
 */
double              t8_forest_element_face_area (t8_forest_t forest,
                                                 t8_locidx_t ltreeid,
                                                 const t8_element_t *element,
                                                 int face);

/** Compute the vertex coordinates of the centroid of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the centroid of \a face.
 * \a forest must be committed when calling this function.
 */
void                t8_forest_element_face_centroid (t8_forest_t forest,
                                                     t8_locidx_t ltreeid,
                                                     const t8_element_t
                                                     *element, int face,
                                                     double centroid[3]);

/** Compute the normal vector of an element's face if a geometry
 * for this tree is registered in the forest's cmesh.
 * Currently implemented for 2D elements only.
 * \param [in]      forest     The forest.
 * \param [in]      ltree_id   The forest local id of the tree in which the element is.
 * \param [in]      element    The element.
 * \param [in]      face       A face of \a element.
 * \param [out]     normal     On output the normal vector of \a element at \a face.
 * \a forest must be committed when calling this function.
 */
void                t8_forest_element_face_normal (t8_forest_t forest,
                                                   t8_locidx_t ltreeid,
                                                   const t8_element_t
                                                   *element, int face,
                                                   double normal[3]);



T8_EXTERN_C_END ();

#endif /* !T8_FOREST_GEOMETRICAL_H */

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

/** \file t8_cad_shape_proximity.hxx
 * Intgrates many CAD and CAM functionalities.
 */

#ifndef T8_CAD_SHAPE_PROXIMITY_HXX
#define T8_CAD_SHAPE_PROXIMITY_HXX

#include <t8.h>
#include <t8_forest.h>

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>
#include <Bnd_OBB.hxx>
#include <NCollection_DefineHArray1.hxx>

/* *INDENT-OFF* */

T8_EXTERN_C_BEGIN ();

typedef NCollection_Array1<Bnd_OBB> Bnd_Array1OfBndOBB; /**< Array to handle oriented bounding boxes (Bnd_OBB) */
DEFINE_HARRAY1(Bnd_HArray1OfBndOBB, Bnd_Array1OfBndOBB) /**< OpenCASCADE memory management */

class t8_cad_shape_proximity
{
public:
  /**
   * Constructor of the cad collsion class. Fills the internal shape with the content of a CAD file.
   * \param [in] filename Path to a CAD file in the .brep or .step format.
   */
  t8_cad_shape_proximity (const char *filename, int use_individual_bbs);

  /**
   * Constructor of the cad class. Fills the internal shape with the given shape.
   * \param [in] shape Occ shape geometry.
   */
  t8_cad_shape_proximity (const TopoDS_Shape shape, int use_individual_bbs);

  /**
   * The destructor
   */
  ~ t8_cad_shape_proximity () {}

  /**
   * Read a brep file and fill internal shape with it.
   * \param [in] fileprefix Prefix of a .brep file from which to extract an occ geometry.
   */
  void
  t8_cad_init (const char *fileprefix, int use_individual_bbs);

  /**
   * Fill the internal shape with the given shape.
   * \param [in] shape Occ shape geometry.
   */
  void
  t8_cad_init (TopoDS_Shape shape, int use_individual_bbs);

  /**
   * Checks if an element is inside the occ shape. Only viable with 
   * axis-oriented hex elements.
   * \param [in] forest          The forest.
   * \param [in] ltreeid         The local tree id of the element.
   * \param [in] element         The element.
   * \return                     0: Element is fully outside of the shape.
   *                             1: Element is partially inside the shape.
   *                             2: Element is fully inside the shape.
   */
  int
  t8_cad_is_element_inside_shape (t8_forest_t forest,
                                  t8_locidx_t ltreeid,
                                  const t8_element_t *element,
                                  int optimize);

  /**
   * Checks if a point is inside the occ shape.
   * \param [in] coords          The coordinates of the point.
   * \param [in] tol             The tolerance.
   * \return                     True if point is inside the shape.
   */
  int
  t8_cad_is_point_inside_shape (const double *coords, double tol, int optimize) const;

protected:
  /**
   * Initializes the internal data structures.
   */
  void                  t8_cad_init_internal_data (const int use_individual_bbs);

  TopoDS_Shape          occ_shape;                                    /**< The OpenCASCADE shape */
  Bnd_OBB               occ_shape_bounding_box;                       /**< Bounding box of the shape */
  Bnd_HArray1OfBndOBB  *occ_shape_individual_bounding_boxes;          /**< Bounding boxes of the subshapes */ 
};

T8_EXTERN_C_END ();

/* *INDENT-ON* */

#endif /* T8_WITH_OCC */

#endif /* !T8_CAD_shape_proximity_HXX */

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

/** \file t8_cad_utils.hxx
 * Utilities for working with CAD data.
 */

#ifndef T8_CAD_UTILS_HXX
#define T8_CAD_UTILS_HXX

#include <t8.h>
#include <t8_eclass.h>

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>

/* *INDENT-OFF* */

/**
 * Reads a CAD file in the BREP, STEP or IGES format and extracts the shape.
 * \param [in] filename  Path to a CAD file in the BREP, STEP or IGES format.
 * \return               The shape of the CAD file. 
 */
TopoDS_Shape
t8_cad_read_cad_file (const char *filename);

/**
 * Uses the diagonal vertices of a hex element to build a
 * cad shape of the element. Faster than \a t8_cad_make_element_shape.
 * \param [in] vertex1   Corner values of an axis-aligned hex element.
 *                       (3 * double)
 * \param [in] vertex2   Values if a corner diagonal of vertex1.
 *                       (3 * double)
 * \return               The cad shape.
 */
TopoDS_Shape
t8_cad_make_axis_aligned_hex_element_shape (const double *vertex1,
                                            const double *vertex2);

/**
 * Uses the vertices of an element to build a cad shape of the ekement.
 * \param [in] vertices  Vertex coordinates of the element in zorder.
 * \param [in] eclass    eclass of the element.
 * \return
 */
TopoDS_Shape
t8_cad_make_element_shape (const double *vertices, const t8_eclass_t eclass);

/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
#endif /* !T8_CAD_UTILS_HXX */

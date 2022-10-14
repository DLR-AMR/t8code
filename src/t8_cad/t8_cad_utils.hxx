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

#if T8_WITH_OCC
#include <TopoDS_Shape.hxx>

/* *INDENT-OFF* */

/**
 * Reads a CAD file in the BRep, STEP or IGES format and extracts the shape.
 * \param [in] filename  Path to a CAD file in the BRep, STEP or IGES format.
 * \return               The shape of the CAD file. 
 */
TopoDS_Shape
t8_cad_read_cad_file (const char *filename);

/* *INDENT-ON* */

#endif /* T8_WITH_OCC */
#endif /* !T8_CAD_UTILS_HXX */

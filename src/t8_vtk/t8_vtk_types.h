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

#ifndef T8_VTK_TYPES
#define T8_VTK_TYPES

#include <t8_eclass.h>

/**
 * Translator between vtk-type of elements and t8code-elements. 
 * Not all elements are supported. Return T8_ECLASS_INVALID for unsupported
 * elements.  
 */
const t8_eclass_t t8_cmesh_vtk_type_to_t8_type[82]
  = { T8_ECLASS_INVALID, T8_ECLASS_VERTEX,  T8_ECLASS_INVALID, T8_ECLASS_LINE,    T8_ECLASS_INVALID, T8_ECLASS_TRIANGLE,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_QUAD,    T8_ECLASS_QUAD,    T8_ECLASS_TET,     T8_ECLASS_HEX,
      T8_ECLASS_HEX,     T8_ECLASS_PRISM,   T8_ECLASS_PYRAMID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID,
      T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID, T8_ECLASS_INVALID };

/**
 * Enumerator for all types of files readable by t8code. 
 */
typedef enum vtk_file_type {
  VTK_FILE_ERROR = -1, /* For Testing purpose. */

  VTK_SERIAL_FILE = 8,
  VTK_UNSTRUCTURED_FILE = VTK_SERIAL_FILE,
  VTK_POLYDATA_FILE,

  VTK_PARALLEL_FILE = 16,
  VTK_PARALLEL_UNSTRUCTURED_FILE = VTK_PARALLEL_FILE, /* For parallel unstructured files. */
  VTK_PARALLEL_POLYDATA_FILE,

  VTK_NUM_TYPES = 5
} vtk_file_type_t;

typedef enum vtk_read_success { read_failure = 0, read_success = 1 } vtk_read_success_t;

#endif /* T8_VTK_TYPES */

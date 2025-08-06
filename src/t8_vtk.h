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

/** file t8_vtk.h
 * This header file collects macros that are needed for
 * the forest and cmesh vtk routines.
 * \see t8_forest_vtk.h \see t8_cmesh_vtk_writer.h \see t8_cmesh_vtk_reader.hxx
 */

#ifndef T8_VTK_H
#define T8_VTK_H

#include <t8.h>

/* typedef and macros */

#define T8_VTK_LOCIDX "Int32"
/* TODO: Paraview has troubles with Int64, so we switch to Int32 and be careful.
 *       Investigate this further. See also vtk macro VTK_USE_64BIT_IDS */
#define T8_VTK_GLOIDX "Int32"

/* TODO: these macros need to be set by configure */
#ifndef T8_VTK_DOUBLES
#define T8_VTK_FLOAT_NAME "Float32"
#define T8_VTK_FLOAT_TYPE float
#else
#define T8_VTK_FLOAT_NAME "Float64"
#define T8_VTK_FLOAT_TYPE double
#endif

#define T8_VTK_FORMAT_STRING "ascii"

#if T8_ENABLE_VTK
#define t8_vtk_locidx_array_type_t vtkTypeInt32Array
#define t8_vtk_gloidx_array_type_t vtkTypeInt64Array
#endif

/* TODO: Add support for integer data type. */
typedef enum {
  T8_VTK_SCALAR, /* One double value per element */
  T8_VTK_VECTOR  /* 3 double values per element */
} t8_vtk_data_type_t;

typedef struct
{
  t8_vtk_data_type_t type;  /**< Describes of which type the data array is */
  char description[BUFSIZ]; /**< String that describes the data. */
  double *data;
  /**< An array of length n*num_local_elements doubles with
                      n = 1 if type = T8_VTK_SCALAR, n = 3 if type = T8_VTK_VECTOR */
} t8_vtk_data_field_t;

T8_EXTERN_C_BEGIN ();

/* function declarations */
/* Writes the pvtu header file that links to the processor local files.
 * It is used by the cmesh and forest vtk routines.
 * This function should only be called by one process.
 * Return T8_SUBROUTINE_SUCCESS on success and T8_SUBROUTINE_FAILURE on failure. */
/* TODO: document */
int
t8_write_pvtu (const char *filename, int num_procs, int write_tree, int write_rank, int write_level, int write_id,
               int num_data, t8_vtk_data_field_t *data);

T8_EXTERN_C_END ();

#endif /* !T8_VTK_H */

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

/** file t8_latlon_data.h
 */

#ifndef T8_LATLON_DATA_H
#define T8_LATLON_DATA_H

#include <t8.h>
#include <t8_forest.h>

/* If we associate data on an X x Y subgrid, it can be 
 * sorted in memory to match the gridcells in different ways.
 */
typedef enum
{
  T8_LATLON_DATA_MESSY,       /* Row-wise storage data[y * X + x] gives data of (x,y). */
  T8_LATLON_DATA_MORTON         /* Morton SFC storage. The data is sorted according to the 
                                 * Morton SFC index in the surrounding quad forest (not the subgrid). */
} T8_LATLON_DATA_NUMBERING;

/* Describes an X by Y subchunk of the grid with
 * data on it. */
typedef struct
{
  const char         *description;  /* The name of this dataset. */
  t8_locidx_t         x_start;      /* Starting x coordinate. */
  t8_locidx_t         y_start;      /* Starting y coordinate. */
  t8_locidx_t         x_length;     /* Number of subgrid cells in x tracer. */
  t8_locidx_t         y_length;     /* Number of subgrid cells in y tracer. */
  t8_locidx_t         z_length;     /* Number of subgrid cells in z tracer. */
  int                 axis;         /* Internal flag to distinguish between different axis configurations e.g. XYZ, YZX, ... */
  int                 x_axis;       /* X axis index in data vector */
  int                 y_axis;       /* Y axis index in data vector */
  int                 z_axis;       /* Z axis index in data vector */
  int                 num_tracers;    /* Dimensionality of the data (1, 2, 3). */
  int                 level;        /* The smallest uniform refinement level of a forest that can have the grid (not the subgrid) as submesh. */
  int                *shape;
  int                tracer_names_size;
  char               *tracer_names;
  double             *data;         /* x_length x y_length x z_length x tracer many data items. For each data item tracer many entries. */
  t8_linearidx_t     *data_ids;     /* Morton index for each grid cells. 
                                     * At first we have (x_lenght x y_length) elements, but when coarsening the number of elements  reduce. */

  double             *data_adapt;
  t8_linearidx_t     *data_ids_adapt;
  T8_LATLON_DATA_NUMBERING numbering; /* Numbering scheme */
} t8_latlon_data_chunk_t;

T8_EXTERN_C_BEGIN ();

/* function declarations */

/**
 * TODO: add doc
 */
t8_latlon_data_chunk_t *
t8_latlon_new_chunk (const char *description, t8_locidx_t x_start, t8_locidx_t y_start,
                     t8_locidx_t x_length, t8_locidx_t y_length, t8_locidx_t z_length,
                     int* shape, int num_tracers, int x_axis, int y_axis, int z_axis, int level,
                     T8_LATLON_DATA_NUMBERING numbering);

/**
 * TODO: add doc
 */
void t8_latlon_chunk_destroy (t8_latlon_data_chunk_t ** pchunk);

/** Given x and y coordinates in an X by Y grid compute
 * the Morton linear id according to a given level of the 
 * element associated with (x,y).
 * 
 * \param [in] x  X-coordinate
 * \param [in] y  Y-coordinate
 * \param [in] level The given refinement level
 * \return        Linear id of the quad the corresponds to the x,y coordinates at \a level.
 */
t8_linearidx_t      t8_latlon_to_linear_id (t8_gloidx_t x, t8_gloidx_t y,
                                            int level);

/** The inverse operation to t8_latlon_to_linear_id.
 * Given a linear id and a refinement level compute the 
 * x and y coordinates of the corresponding grid cell.
 * \param [in] linear_id  The Morton linear id of a quad element.
 * \param [in] level      The considered refinement level.
 * \param [out] x         On output filled with the grid x coordinate.
 * \param [out] y         On output filled with the grid y coordinate.
 */
void                t8_latlon_linear_id_to_latlon (t8_linearidx_t linear_id,
                                                   int level, t8_gloidx_t * x,
                                                   t8_gloidx_t * y);

int t8_latlon_get_tracer_idx(t8_latlon_data_chunk_t * data_chunk, char* tracer, bool add_if_missing);

void  t8_latlon_set_dimension(t8_latlon_data_chunk_t * data_chunk, char* tracer, double**** data);

void
t8_latlon_data_apply_morton_order (t8_forest_t *forest, t8_latlon_data_chunk_t * data_chunk);


/* Create a data chunk with given num_tracers and numbering,
 * fill it with data and then change the numbering to Morton.
 */
void                t8_latlon_data_test (t8_locidx_t x_start,
                                         t8_locidx_t y_start,
                                         t8_locidx_t x_length,
                                         t8_locidx_t y_length,
                                         int* shape,
                                         int tracer, int x_axis, 
                                         int y_axis, int z_axis, int level,
                                         T8_LATLON_DATA_NUMBERING
                                         numbering,
                                         t8_gloidx_t x_length_global,
                                         t8_gloidx_t y_length_global);

T8_EXTERN_C_END ();

#endif /* !T8_LATLON_DATA_H */

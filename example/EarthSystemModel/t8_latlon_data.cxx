/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <sc_options.h>
#include <sc_refcount.h>
#include <t8.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common_cxx.hxx>
#include <t8_forest.h>
#include <t8_cmesh_vtk.h>
#include <t8_schemes/t8_default/t8_default_quad_cxx.hxx>
#include <p4est_bits.h>
#include "t8_latlon_data.h"

/* Given x and y coordinates in a X by Y grid compute the quad
 * element that contains these coordinates when the grid is embedded
 * in the lower left of a quad tree.
 *  _ _ _ _ _ _
 * |           |
 * |           |  X = 3, Y = 2  
 * |_____ __   |  * = (0, 1) x and y coordinates
 * |_*|__|__|  |
 * |__|__|__|__|
 * 
 */
static void
t8_latlon_to_element (t8_gloidx_t x, t8_gloidx_t y, int level,
                      t8_element_t * quad_element)
{
  p4est_quadrant_t   *quad;
  t8_default_scheme_quad_c quad_scheme;

  quad = (p4est_quadrant_t *) (quad_element);
  quad->level = level;
  quad->x = x << (P4EST_MAXLEVEL - level);
  quad->y = y << (P4EST_MAXLEVEL - level);
}

t8_linearidx_t
t8_latlon_to_linear_id (t8_gloidx_t x, t8_gloidx_t y, int level)
{
  t8_element_t       *elem;
  t8_linearidx_t      linear_id;
  t8_default_scheme_quad_c quad_scheme;

  quad_scheme.t8_element_new (1, &elem);

  t8_latlon_to_element (x, y, level, elem);
  linear_id = quad_scheme.t8_element_get_linear_id (elem, level);
  quad_scheme.t8_element_destroy (1, &elem);

  return linear_id;
}

void
t8_latlon_linear_id_to_latlon (t8_linearidx_t linear_id, int level,
                               t8_gloidx_t * x, t8_gloidx_t * y)
{
  t8_element_t       *elem;
  p4est_quadrant_t   *quad;
  t8_default_scheme_quad_c quad_scheme;
  quad_scheme.t8_element_new (1, &elem);

  quad_scheme.t8_element_set_linear_id (elem, level, linear_id);

  quad = (p4est_quadrant_t *) (elem);
  *x = quad->x >> (P4EST_MAXLEVEL - level);
  *y = quad->y >> (P4EST_MAXLEVEL - level);

  quad_scheme.t8_element_destroy (1, &elem);
}

void
t8_latlon_data_index_to_latloan (t8_latlon_data_chunk_t * data_chunk,
                                 t8_locidx_t array_index,
                                 t8_gloidx_t * x_coord, t8_gloidx_t * y_coord)
{
  t8_locidx_t         x_offset;
  t8_locidx_t         y_offset;
  t8_linearidx_t      data_id;
  switch (data_chunk->numbering) {
  case T8_LATLON_DATA_XSTRIPE:
    /*   _ _ _ _
     *  |   ____|
     *  |  |--->| 
     *  |  |--->|
     *  |_ _ _ _|
     */
    /* Compute the x and y position within the data chunk. */
    x_offset = array_index % data_chunk->x_length;
    y_offset = array_index / data_chunk->x_length;
    /* Add the start indices of the chunk to compute the x and y
     * position in the grid. */
    *x_coord = data_chunk->x_start + x_offset;
    *y_coord = data_chunk->y_start + y_offset;
    return;
  case T8_LATLON_DATA_YSTRIPE:
    /*   _ _ _ _
     *  |   ____|
     *  |  |^^^^| 
     *  |  ||||||
     *  |_ _ _ _|
     */
    /* Compute the x and y position within the data chunk. */
    x_offset = array_index / data_chunk->y_length;
    y_offset = array_index % data_chunk->y_length;
    /* Add the start indices of the chunk to compute the x and y
     * position in the grid. */
    *x_coord = data_chunk->x_start + x_offset;
    *y_coord = data_chunk->y_start + y_offset;
    return;
  case T8_LATLON_DATA_MORTON:
    data_id = data_chunk->data_ids[array_index];
    t8_latlon_linear_id_to_latlon (data_id, data_chunk->level, x_coord,
                                   y_coord);
    return;
  default:
    /* Safety measure to prevent us from adding numbering schemes
     * in the future and not adding them here. */
    SC_ABORT_NOT_REACHED ();
  }
}

t8_latlon_data_chunk_t *
t8_latlon_new_chunk (t8_locidx_t x_start, t8_locidx_t y_start,
                     t8_locidx_t x_length, t8_locidx_t y_length,
                     int dimension, int level,
                     T8_LATLON_DATA_NUMBERING numbering,
                     const char *description)
{
  t8_latlon_data_chunk_t *chunk = T8_ALLOC (t8_latlon_data_chunk_t, 1);

  /* Since we currently do not know how to compute the morton ids
   * beforehand, we cannot fill a chunk in morton order without knowing
   * also the ids. */
  T8_ASSERT (numbering != T8_LATLON_DATA_MORTON);
  chunk->x_start = x_start;
  chunk->y_start = y_start;
  chunk->x_length = x_length;
  chunk->y_length = y_length;
  chunk->numbering = numbering;
  chunk->dimension = dimension;
  chunk->data = T8_ALLOC (double, dimension * x_length * y_length);
  chunk->data_ids = NULL;
  chunk->description = description;
  chunk->level = level;

  return chunk;
}

void
t8_latlon_chunk_destroy (t8_latlon_data_chunk_t ** pchunk)
{
  T8_ASSERT (pchunk != NULL);
  t8_latlon_data_chunk_t *chunk = *pchunk;
  T8_ASSERT (chunk != NULL);

  /* Free all data and the chunk itself. */
  T8_FREE (chunk->data);
  T8_FREE (chunk->data_ids);
  T8_FREE (chunk);

  /* Set the input pointer to NULL */
  *pchunk = NULL;
}

/* Given two indices into an array data_ids return
 * < 0 if data_ids[index1] < data_ids[index2]
 *   0 if data_ids[index1] == data_ids[index2]
 * > 0 if data_ids[index1] > data_ids[index2]
 * 
 * We use this to sort an array of indices according to the values in
 * data_ids. */
static int
t8_latlon_compare_indices (const void *index1, const void *index2,
                           void *data_ids)
{
  t8_linearidx_t     *data_ids_idx = (t8_linearidx_t *) data_ids;
  size_t              index1_value = *(size_t *) index1;
  size_t              index2_value = *(size_t *) index2;

  /* We cannot use the 'return a - b' trick since t8_linearidx_t may be 
   * of unsigned integer type (at time of writing this it is uint64_t) */
  return data_ids_idx[index1_value] < data_ids_idx[index2_value]
    ? -1 : data_ids_idx[index1_value] == data_ids_idx[index2_value]
    ? 0 : 1;
}

/* Change the numebring of a data chunk to morton numbering. */
void
t8_latlon_data_apply_morton_order (t8_latlon_data_chunk_t * data_chunk)
{
  t8_locidx_t         num_grid_elements =
    data_chunk->x_length * data_chunk->y_length;
  if (data_chunk->numbering == T8_LATLON_DATA_MORTON) {
    /* This data is already in Morton order. */
    return;
  }
  size_t             *permutation;
  /* Allocate array to store morton indices of the data items. */
  T8_ASSERT (data_chunk->data_ids == NULL);
  data_chunk->data_ids = T8_ALLOC (t8_linearidx_t, num_grid_elements);
  permutation = T8_ALLOC (size_t, num_grid_elements);
  /* Compute linear ids for all indices. 
   * Fill permutation array with 0, 1, 2, 3, ... */
  {
    t8_locidx_t         index;
    t8_gloidx_t         x_coord, y_coord;
    t8_debugf ("Building Morton ids:\n");
    for (index = 0; index < num_grid_elements; ++index) {
      /* Compute x and y coordinate in whole grid of this array index. */
      t8_latlon_data_index_to_latloan (data_chunk, index, &x_coord, &y_coord);
      data_chunk->data_ids[index] =
        t8_latlon_to_linear_id (x_coord, y_coord, data_chunk->level);
      permutation[index] = index;
      t8_debugf (" %li\n", data_chunk->data_ids[index]);
    }
  }
  /* We now sort the data according to the data_ids.
   * We do this by sorting the permutation array and them applying
   * this permutation to the data_ids and data arrays. */
  qsort_r (permutation, num_grid_elements, sizeof (size_t),
           t8_latlon_compare_indices, data_chunk->data_ids);
  /* Now the correct order of the data_ids is
   * data_ids[permutation[0]]
   * data_ids[permutation[1]]
   * data_ids[permutation[2]]
   *    ...
   */
  {
    t8_locidx_t         index;
    t8_debugf ("newindices:\n");
    for (index = 0; index < num_grid_elements; ++index) {
      t8_debugf ("%zd\n", permutation[index]);
    }
  }
  {
    /* Reorder the data and the ids.
     * To do so, we create a new array and copy the data 
     * over in the correct order. */
    /* NOTE: We could use sc_array_permute, but this would require the inverse of
     *       permutation and currently we did not figure out how to compute it.
     */
    t8_locidx_t         index;
    int                 dim = data_chunk->dimension;
    double             *data_new = T8_ALLOC (double, dim * num_grid_elements);
    t8_linearidx_t     *data_ids_new =
      T8_ALLOC (t8_linearidx_t, num_grid_elements);
    /* Copy the data */
    for (index = 0; index < num_grid_elements; ++index) {
      data_ids_new[index] = data_chunk->data_ids[permutation[index]];
      /* Copy dim many doubles over */
      memcpy (data_new + dim * index,
              data_chunk->data + dim * permutation[index],
              dim * sizeof (double));
    }
    /* Replace the original arrays */
    T8_FREE (data_chunk->data);
    T8_FREE (data_chunk->data_ids);
    data_chunk->data_ids = data_ids_new;
    data_chunk->data = data_new;
  }

  T8_FREE (permutation);
}

/* Given a subgrid in a global grid that is represented in a partitioned forest,
 * determine those processes that have elements in the forest that lie in the
 * subgrid.
 * The subgrid data will be ordered in Morton order (if not already on input)
 * and we will return the indices of the process bounds.
 *         _ _ _ _ _ _
 * data:  |_|_|_|_|_|_|
 * 
 * procs:  p_0|p_1  | p_3
 * 
 * Will create 
 *  num_processes: 3
 *  process_bounds: 0, 2, 5, 6
 *  process_ranks:  0, 1, 3
 */
void
t8_latlon_data_determine_process_bounds (t8_forest_t forest,
                                         t8_latlon_data_chunk_t * chunk,
                                         t8_gloidx_t x_length_global,
                                         t8_gloidx_t y_length_global,
                                         int *num_processes,
                                         size_t * process_offsets,
                                         int *process_ranks)
{
  if (chunk->numbering != T8_LATLON_DATA_MORTON) {
    /* Apply Morton order if not already set. */
    t8_latlon_data_apply_morton_order (chunk);
  }
  /* TODO: We need to perform a multiple binary search to identify
   * the process borders.
   * Alternatively, we could try to find them while we establish
   * the Morton order and copy the arrays in the correct order.
   */
}

void
t8_latlon_data_test (t8_locidx_t x_start, t8_locidx_t y_start,
                     t8_locidx_t x_length, t8_locidx_t y_length,
                     int dimension, int level,
                     T8_LATLON_DATA_NUMBERING numbering,
                     t8_gloidx_t x_length_global, t8_gloidx_t y_length_global)
{
  t8_latlon_data_chunk_t *chunk =
    t8_latlon_new_chunk (x_start, y_start, x_length, y_length, dimension,
                         level, numbering, "test_chunk");
  t8_locidx_t         num_grid_items = x_length * y_length;
  T8_ASSERT (numbering != T8_LATLON_DATA_MORTON);
  {
    /*  Fill with values 0, 1, 2 , ...  */
    t8_locidx_t         index;
    for (index = 0; index < num_grid_items; ++index) {
      chunk->data[dimension * index] = dimension * index;
      if (dimension > 1) {
        chunk->data[dimension * index + 1] = dimension * index;
      }
      if (dimension > 2) {
        chunk->data[dimension * index + 2] = dimension * index;
      }
    }
  }
  {
    /* Print data. */
    t8_debugf ("Dataset %s with (%i,%i,%i,%i) and numbering %i:\n",
               chunk->description, chunk->x_start, chunk->y_start,
               chunk->x_length, chunk->y_length, chunk->numbering);
    t8_locidx_t         index;
    for (index = 0; index < num_grid_items; ++index) {
      t8_debugf ("%.2f\n", chunk->data[index]);
    }
    t8_latlon_data_apply_morton_order (chunk);
    t8_debugf ("Applied Morton order:\n");
    for (index = 0; index < num_grid_items; ++index) {
      t8_debugf ("%.2f %.2f %.2f (%li)\n", chunk->data[dimension * index],
                 dimension > 1 ? chunk->data[dimension * index + 1] : -1,
                 dimension > 2 ? chunk->data[dimension * index + 2] : -1,
                 chunk->data_ids[index]);
    }
  }
  t8_latlon_chunk_destroy (&chunk);
}

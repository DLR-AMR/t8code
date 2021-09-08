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


#define XYZ 6
#define XZY 9
#define YXZ 18
#define YZX 24
#define ZXY 33
#define ZYX 36


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
 * 
 * Note: This is essentially the invers of the t8_latlon_refine_grid_cuts_elements check.
 * 
 * @author Holke Johannes
 */
static void
t8_latlon_to_element (t8_gloidx_t x, t8_gloidx_t y, int level,
                      t8_element_t * quad_element)
{
  p4est_quadrant_t   *quad;
  t8_default_scheme_quad_c quad_scheme;

  quad = (p4est_quadrant_t *) (quad_element);
  /* The level of the quad is the given level.
   * Its x and y coordinates are the given x and y coordinates in the grid,
   * shifted to be in the range [0,2^L], L = P4EST_MAXLEVEL. */
  quad->level = level;
  quad->x = x << (P4EST_MAXLEVEL - level);
  quad->y = y << (P4EST_MAXLEVEL - level);
}

/* Given x and y coordinates in an X by Y grid compute
 * the Morton linear id according to a given level of the 
 * element associated with (x,y).
 */
t8_linearidx_t
t8_latlon_to_linear_id (t8_gloidx_t x, t8_gloidx_t y, int level)
{
  t8_element_t       *elem;
  t8_linearidx_t      linear_id;
  t8_default_scheme_quad_c quad_scheme;

  /* Allocate memory for a new quad */
  quad_scheme.t8_element_new (1, &elem);

  /* Initialize the quad to match the x,y coordinates */
  t8_latlon_to_element (x, y, level, elem);
  /* Compute the linear id of the quad. */
  linear_id = quad_scheme.t8_element_get_linear_id (elem, level);

  /* Destroy the quad. */
  quad_scheme.t8_element_destroy (1, &elem);

  /* Return the linear id. */
  return linear_id;
}

/* The inverse operation to t8_latlon_to_linear_id.
 * Given a linear id and a refinement level compute the 
 * x and y coordinates of the corresponding grid cell. 
 * 
 * @author Holke Johannes
 */
void
t8_latlon_linear_id_to_latlon (t8_linearidx_t linear_id, int level,
                               t8_gloidx_t * x, t8_gloidx_t * y)
{
  t8_element_t       *elem;
  p4est_quadrant_t   *quad;
  t8_default_scheme_quad_c quad_scheme;

  /* Allocate memory for a new quad */
  quad_scheme.t8_element_new (1, &elem);

  /* Initialize it according to the linear id. */
  quad_scheme.t8_element_set_linear_id (elem, level, linear_id);

  /* Compute x and y from the linear id. */
  quad = (p4est_quadrant_t *) (elem);
  *x = quad->x >> (P4EST_MAXLEVEL - level);
  *y = quad->y >> (P4EST_MAXLEVEL - level);

  /* Clean-up memory. */
  quad_scheme.t8_element_destroy (1, &elem);
}

/* Given a data chunk (in MESSy or Morton order)
 * and an array index into its data array compute the corresponding
 * x and y coordinates in the grid. 
 * 
 * @author Holke Johannes, Spataro Luca
 */
void
t8_latlon_data_index_to_latlon (t8_latlon_data_chunk_t * data_chunk,
                                 t8_locidx_t array_index,
                                 t8_gloidx_t * x_coord, t8_gloidx_t * y_coord)
{
  t8_locidx_t         x_offset;
  t8_locidx_t         y_offset;
  t8_linearidx_t      data_id;
  switch (data_chunk->numbering) {
  case T8_LATLON_DATA_MESSY:
    /* row strides */
    /* Compute the x and y position within the data chunk. */
    x_offset = array_index % data_chunk->x_length;
    y_offset = array_index / data_chunk->x_length;
    
    /* Add the start indices of the chunk to compute the x and y
    * position in the grid. */
    *x_coord = data_chunk->x_start + x_offset;
    *y_coord = data_chunk->y_start + y_offset;
    return;
  case T8_LATLON_DATA_MORTON:
    /* Get the Morton index of the corresponding element. */
    data_id = data_chunk->data_ids[array_index];
    /* Compute the lat/lon coordinates from the index. */
    t8_latlon_linear_id_to_latlon (data_id, data_chunk->level, x_coord,
                                   y_coord);
    return;
  default:
    /* Safety measure to prevent us from adding numbering schemes
     * in the future and not adding them here. */
    SC_ABORT_NOT_REACHED ();
  }
}

/* 
 * Create a new data chunk with given num_tracers and numbering. 
 * 
 * @author Holke Johannes, Spataro Luca
 */
t8_latlon_data_chunk_t *
t8_latlon_new_chunk (const char *description, t8_locidx_t x_start, t8_locidx_t y_start,
                     t8_locidx_t x_length, t8_locidx_t y_length, t8_locidx_t z_length,
                     int *shape, int num_tracers, int x_axis, int y_axis, int z_axis, int level,
                     double missing_value, T8_LATLON_DATA_NUMBERING numbering)
{
  t8_latlon_data_chunk_t *chunk = T8_ALLOC (t8_latlon_data_chunk_t, 1);

  /* Since we currently do not know how to compute the morton ids
   * beforehand, we cannot fill a chunk in morton order without knowing
   * also the ids. */
  T8_ASSERT (numbering != T8_LATLON_DATA_MORTON);

  /* Initialize all values and allocate data array */
  chunk->description = description;
  chunk->level = level;
  chunk->x_start = x_start;
  chunk->y_start = y_start;
  chunk->x_length = x_length;
  chunk->y_length = y_length;
  chunk->z_length = z_length;
  chunk->shape = shape;
  chunk->numbering = numbering;
  chunk->num_tracers = num_tracers;
  chunk->x_axis = x_axis;
  chunk->y_axis = y_axis;
  chunk->z_axis = z_axis;
  chunk->missing_value = missing_value;

  chunk->data = T8_ALLOC_ZERO(double, x_length * y_length * z_length * num_tracers);
  chunk->data_ids = T8_ALLOC_ZERO(t8_linearidx_t, x_length * y_length);
  chunk->data_adapt = NULL;
  chunk->data_ids_adapt = NULL;
  
  chunk->tracer_names_size = 0;
  chunk->tracer_names = T8_ALLOC(char, BUFSIZ * num_tracers);
  
  /* bit concatenate axis configuration
   * e.g. x = 0 = 00, y = 1 = 01, z = 2 = 10 => 000110
   *      x = 2 = 10, y = 0 = 00, z = 1 = 01 => 100001
   * ...
   */
  chunk->axis = x_axis << 4 | y_axis << 2 | z_axis;

  return chunk;
}

/* Destroy a data chunk 
 *
 * @author Holke Johannes, Spataro Luca
 */
void
t8_latlon_chunk_destroy (t8_latlon_data_chunk_t ** pchunk)
{
  T8_ASSERT (pchunk != NULL);
  t8_latlon_data_chunk_t *chunk = *pchunk;
  T8_ASSERT (chunk != NULL);

  /* Free all data and the chunk itself. */
  T8_FREE (chunk->data);
  T8_FREE (chunk->data_ids);

  T8_FREE (chunk->tracer_names);

  T8_FREE (chunk->shape);

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
 * data_ids. 
 * 
 * @author Holke Johannes
 * */
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

/* Retrive value from input data on given (x_coord, y_coord, z_coord, dim) coordinate 

 * @author Spataro Luca
 */
double t8_latlon_get_dimension_value(int axis, double ****data, int x_coord, 
                             int y_coord, int z_coord, int dimension) {
  double value;
  switch(axis) {
    case XYZ:
    value = data[x_coord][y_coord][z_coord][dimension];
    break;
    case XZY:
    value = data[x_coord][z_coord][y_coord][dimension];
    break;
    case YXZ:
    value = data[y_coord][x_coord][z_coord][dimension];
    break;
    case YZX:
    value = data[y_coord][z_coord][x_coord][dimension];
    break;
    case ZXY:
    value = data[z_coord][x_coord][y_coord][dimension];
    break;
    case ZYX:
    value = data[z_coord][y_coord][x_coord][dimension];
    break;
    /* TODO: Add default label with SC_ABORT */
  }

  return value;
}

int t8_latlon_get_tracer_idx(t8_latlon_data_chunk_t * data_chunk, char* dimension, bool add_if_missing) {
  int idx;
  
  /* search for dimension name */
  for(idx = 0; idx < data_chunk->tracer_names_size; ++idx) {
    if (strncmp(dimension, data_chunk->tracer_names + idx * BUFSIZ, BUFSIZ) == 0) {
      return idx;
    }
  }

  if(add_if_missing) {
    /* not found, so add to end */
    t8_debugf("dimension %s missing, adding it to list \n", dimension);
    idx = data_chunk->tracer_names_size;
    strncpy(data_chunk->tracer_names + idx * BUFSIZ, dimension, BUFSIZ);
    data_chunk->tracer_names_size = (data_chunk->tracer_names_size) + 1;
    return idx;
  }
  
  return -1;
}

void  t8_latlon_set_dimension(t8_latlon_data_chunk_t * data_chunk, char* tracer, double**** data) {

}

/* set dimension value inside input data on given x_coord, y_coord coordinate 

 * @author Spataro Luca
 */
void t8_latlon_set_dimension_value(int axis, double ****data, int x_coord, 
                             int y_coord, int z_coord, int dimension, double value) {
  switch(axis) {
    case XYZ:
    data[x_coord][y_coord][z_coord][dimension] = value;
    break;
    case XZY:
    data[x_coord][z_coord][y_coord][dimension] = value;
    break;
    case YXZ:
    data[y_coord][x_coord][z_coord][dimension] = value;
    break;
    case YZX:
    data[y_coord][z_coord][x_coord][dimension] = value;
    break;
    case ZXY:
    data[z_coord][x_coord][y_coord][dimension] = value;
    break;
    case ZYX:
    data[z_coord][y_coord][x_coord][dimension] = value;
    break;
  }
}

/* Change the numbering of a data chunk to morton numbering. 
 * 
 * @author Spataro Luca
 */
void
t8_latlon_data_apply_morton_order (t8_forest_t *forest, t8_latlon_data_chunk_t * data_chunk)
{
  t8_debugf ("Applying morton order\n");
  
  if (data_chunk->numbering == T8_LATLON_DATA_MORTON) {
    /* This data is already in Morton order. */
    return;
  }

  int z_length = data_chunk->z_length;
  int element_length = z_length * data_chunk->num_tracers;
  int num_elements = t8_forest_get_local_num_elements(*forest);
  int num_data_elements = num_elements * element_length;

  int coords[2];
  int old;
  t8_gloidx_t x, y, index;

  t8_element_t *elem;
  t8_default_scheme_quad_c quad_scheme;

  t8_debugf("num_data_elements: %d, num_element: %d \n", num_data_elements, num_elements);

  double             *data_new = T8_ALLOC (double, num_data_elements);
  t8_linearidx_t     *data_ids_new = T8_ALLOC_ZERO (t8_linearidx_t, num_elements);

  for (index = 0; index < num_elements; ++index) {
    elem = t8_forest_get_element_in_tree(*forest, 0, index);
    quad_scheme.t8_element_vertex_coords(elem, 0, coords);
    
    x = coords[0] >> (P4EST_MAXLEVEL - data_chunk->level);
    y = coords[1] >> (P4EST_MAXLEVEL - data_chunk->level);
    if(x < data_chunk->x_length && y < data_chunk->y_length) {
      old = y * data_chunk->x_length * element_length + x * element_length;

      /* Copy dim many doubles over */
      memcpy (data_new + index * element_length,
              data_chunk->data + old,
              element_length * sizeof (double));
    } else {
      for(int i=0; i<element_length; ++i) {
        data_new[index * element_length + i] = data_chunk->missing_value;
      }
    }
  }

  /* Replace the original arrays */
  T8_FREE (data_chunk->data);
  T8_FREE (data_chunk->data_ids);
  
  data_chunk->data_ids = data_ids_new;
  data_chunk->data = data_new;

  data_chunk->numbering = T8_LATLON_DATA_MORTON;

  t8_debugf ("Morton order applied\n");
}


/* WIP: Given a subgrid in a global grid that is represented in a partitioned forest,
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
    t8_latlon_data_apply_morton_order (NULL, chunk);
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
                      int* shape,
                     int dimension, int x_axis, int y_axis, int z_axis,
                     int level, T8_LATLON_DATA_NUMBERING numbering,
                     t8_gloidx_t x_length_global, t8_gloidx_t y_length_global)
{
  T8_ASSERT (x_start + x_length <= x_length_global);
  T8_ASSERT (y_start + y_length <= y_length_global);
  t8_latlon_data_chunk_t *chunk =
    t8_latlon_new_chunk ("test_chunk", x_start, y_start, x_length, y_length, 1, shape, dimension,
                         x_axis, y_axis, z_axis, level, 0.0, numbering);
  t8_locidx_t         num_grid_items = x_length * y_length;
  T8_ASSERT (numbering != T8_LATLON_DATA_MORTON);
  {
    int d1 = x_axis == 0 ? x_length : (y_axis == 0 ? y_length : 1);
    int d2 = x_axis == 1 ? x_length : (y_axis == 1 ? y_length : 1);
    int d3 = x_axis == 2 ? x_length : (y_axis == 2 ? y_length : 1);
    /*  Fill with values 0, 1, 2 , ...  */
    t8_locidx_t         index;
    
    /*for (index = 0; index < num_grid_items; ++index) {
      x = index % d1;
      y = index / d2;

      chunk->in[d1][d2][d3] = dimension * index;
      chunk->data[dimension * index] = dimension * index;
      if (dimension > 1) {
        chunk->data[dimension * index + 1] = dimension * index;
      }
      if (dimension > 2) {
        chunk->data[dimension * index + 2] = dimension * index;
      }
    }*/
  }
  {
    /* Print data. */
    t8_debugf ("Dataset %s with (%i,%i,%i,%i) and numbering %i:\n",
               chunk->description, chunk->x_start, chunk->y_start,
               chunk->x_length, chunk->y_length, chunk->numbering);
    t8_locidx_t         index;
    /*for (index = 0; index < num_grid_items; ++index) {
      t8_debugf ("%.2f\n", chunk->data[index]);
    }*/
    /* Change to Morton order */
    t8_latlon_data_apply_morton_order (NULL, chunk);
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

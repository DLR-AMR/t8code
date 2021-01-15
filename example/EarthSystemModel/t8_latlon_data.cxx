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
    /* Check in what format the data is stored */
    if(data_chunk->x_axis < data_chunk->y_axis) {
      /* column strides */
      /* Compute the x and y position within the data chunk. */
      x_offset = array_index / data_chunk->y_length;
      y_offset = array_index % data_chunk->y_length;
    } else {
      /* row strides */
      /* Compute the x and y position within the data chunk. */
      x_offset = array_index % data_chunk->x_length;
      y_offset = array_index / data_chunk->x_length;
    }
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
 * Create a new data chunk with given dimensions and numbering. 
 * 
 * @author Holke Johannes, Spataro Luca
 */
t8_latlon_data_chunk_t *
t8_latlon_new_chunk (const char *description, t8_locidx_t x_start, t8_locidx_t y_start,
                     t8_locidx_t x_length, t8_locidx_t y_length, t8_locidx_t z_length,
                     int dimensions, int x_axis, int y_axis, int z_axis, int level,
                     T8_LATLON_DATA_NUMBERING numbering
                     )
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
  chunk->numbering = numbering;
  chunk->dimensions = dimensions;
  chunk->x_axis = x_axis;
  chunk->y_axis = y_axis;
  chunk->z_axis = z_axis;

  chunk->data = NULL;
  chunk->data_ids = NULL;
  chunk->data_adapt = NULL;
  chunk->data_ids_adapt = NULL;
  
  chunk->dimension_names_size = 0;
  chunk->dimension_names = T8_ALLOC(char*, dimensions);

  for(int d = 0; d < dimensions; ++d) {
    chunk->dimension_names[d] = T8_ALLOC(char, sizeof(char) * BUFSIZ);
  }

  /* bit concatenate axis configuration
   * e.g. x = 0 = 00, y = 1 = 01, z = 2 = 10 => 000110
   *      x = 2 = 10, y = 0 = 00, z = 1 = 01 => 100001
   * ...
   */
  chunk->axis = x_axis << 4 | y_axis << 2 | z_axis;

  /* determine axis length */
  int len_1 = x_axis == 0 ? x_length : (y_axis == 0 ? y_length : z_length);
  int len_2 = x_axis == 1 ? x_length : (y_axis == 1 ? y_length : z_length);
  int len_3 = x_axis == 2 ? x_length : (y_axis == 2 ? y_length : z_length);

  /* allocate input data array */
  chunk->in = T8_ALLOC(double***, len_1);
  int i = 0, j = 0, k = 0;
  for (i = 0; i < len_1; ++i) {
      chunk->in[i] = T8_ALLOC(double**, len_2);
      for (j = 0; j < len_2; ++j) {
          chunk->in[i][j] = T8_ALLOC(double*, len_3);
          for(k = 0; k < len_3; ++k) {
            chunk->in[i][j][k] = T8_ALLOC(double, dimensions);
          }
      }
  }

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

  T8_FREE (chunk->dimension_names);

  int d1 = chunk->x_axis == 0 ? chunk->x_length : (chunk->y_axis == 0 ? chunk->y_length : chunk->z_length);
  int d2 = chunk->x_axis == 1 ? chunk->x_length : (chunk->y_axis == 1 ? chunk->y_length : chunk->z_length);
  int d3 = chunk->x_axis == 2 ? chunk->x_length : (chunk->y_axis == 2 ? chunk->y_length : chunk->z_length);

  /* free input data */
  int i = 0, j = 0, k = 0;
  for (i = 0; i < d1; ++i) {
      for (j = 0; j < d2; ++j) {
          for(k = 0; k < d3; ++k) {
            T8_FREE(chunk->in[i][j][k]);
          }
          T8_FREE(chunk->in[i][j]);
      }
      T8_FREE(chunk->in[i]);
  }
  T8_FREE(chunk->in);

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
  }

  return value;
}

int t8_latlon_get_dimension_idx(t8_latlon_data_chunk_t * data_chunk, char* dimension, bool add_if_missing) {
  int idx;
  /* search for dimension name */
  for(idx = 0; idx < data_chunk->dimension_names_size; ++idx) {
    if (strncmp(dimension, (data_chunk->dimension_names[idx]), BUFSIZ) == 0) {
      return idx;
    }
  }

  if(add_if_missing) {
    /* not found, so add to end */
    idx = data_chunk->dimension_names_size;
    (data_chunk->dimension_names)[idx] = dimension;
    data_chunk->dimension_names_size = (data_chunk->dimension_names_size) + 1;

    return idx;
  }

  return -1;
}

void  t8_latlon_set_dimension(t8_latlon_data_chunk_t * data_chunk, char* dimension_name, double**** data) {
  int axis = data_chunk->axis;
  int x, y, z;
  double value;

  int dimension = t8_latlon_get_dimension_idx(data_chunk, dimension_name, true);

  t8_debugf("add dimension %s at index %d\n", dimension_name, dimension);

  /* maybe we should iterate in the correct axis order */
  for(x=0; x < data_chunk->x_length; ++x) {
    for(y=0; y < data_chunk->y_length; ++y) {
      for(z=0; z < data_chunk->z_length; ++z) {
        value = t8_latlon_get_dimension_value(axis, data, x, y, z, 0);
        t8_latlon_set_dimension_value(axis, data_chunk->in, x, y, z, dimension, value);
      }
    }
  }
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
 * @author Holke Johannes, Spataro Luca
 */
void
t8_latlon_data_apply_morton_order (t8_latlon_data_chunk_t * data_chunk)
{
  t8_debugf ("Applying morton order\n");

  if (data_chunk->numbering == T8_LATLON_DATA_MORTON) {
    /* This data is already in Morton order. */
    return;
  }
  size_t             *permutation;
  /* Allocate array to store morton indices of the data items. */
  T8_ASSERT (data_chunk->data_ids == NULL);
  int num_dimension = data_chunk->dimensions;
  int z_length = data_chunk->z_length;
  int element_length = z_length * num_dimension;
  int num_grid_elements = data_chunk->x_length * data_chunk->y_length;
  int num_data_elements = num_grid_elements * element_length;

  permutation = T8_ALLOC (size_t, num_grid_elements);
  data_chunk->data_ids = T8_ALLOC (t8_linearidx_t, num_grid_elements);
  data_chunk->data = T8_ALLOC (double, num_data_elements);

  /* Compute linear ids for all indices. 
   * Fill permutation array with 0, 1, 2, 3, ... */
  {
    t8_locidx_t         index;
    t8_locidx_t         d, z, i;
    t8_gloidx_t         x_coord, y_coord;
    t8_debugf ("Building Morton ids:\n");
    for (index = 0; index < num_grid_elements; ++index) {
      /* Compute x and y coordinate in whole grid for array index. */
      t8_latlon_data_index_to_latlon (data_chunk, index, &x_coord, &y_coord);

      /* Retrive value from input data */
      for (z = 0; z < z_length; ++z) {
        for (d = 0; d < num_dimension; ++d) {
          i = index * element_length + z * num_dimension + d;
          data_chunk->data[i] = t8_latlon_get_dimension_value(data_chunk->axis, data_chunk->in, x_coord, y_coord, z, d);
          t8_debugf ("(%d)[%d, %d, %d][%d]: %.2f\n", index, x_coord, y_coord, z, d, data_chunk->data[i]);
        }
      }

      data_chunk->data_ids[index] =
        t8_latlon_to_linear_id (x_coord, y_coord, data_chunk->level);
      permutation[index] = index;
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
    double             *data_new = T8_ALLOC (double, num_data_elements);
    t8_linearidx_t     *data_ids_new = T8_ALLOC (t8_linearidx_t, num_grid_elements);
    /* Copy the data */
    for (index = 0; index < num_grid_elements; ++index) {
      data_ids_new[index] = data_chunk->data_ids[permutation[index]];
      /* Copy dim many doubles over */
      memcpy (data_new + index * element_length,
              data_chunk->data + permutation[index] * element_length,
              element_length * sizeof (double));
    }
    /* Replace the original arrays */
    T8_FREE (data_chunk->data);
    T8_FREE (data_chunk->data_ids);
    data_chunk->data_ids = data_ids_new;
    data_chunk->data = data_new;
  }

  T8_FREE (permutation);
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
                     int dimension, int x_axis, int y_axis, int z_axis, 
                     int level, T8_LATLON_DATA_NUMBERING numbering,
                     t8_gloidx_t x_length_global, t8_gloidx_t y_length_global)
{
  T8_ASSERT (x_start + x_length <= x_length_global);
  T8_ASSERT (y_start + y_length <= y_length_global);
  t8_latlon_data_chunk_t *chunk =
    t8_latlon_new_chunk ("test_chunk", x_start, y_start, x_length, y_length, 1, dimension,
                         x_axis, y_axis, z_axis, level, numbering);
  t8_locidx_t         num_grid_items = x_length * y_length;
  T8_ASSERT (numbering != T8_LATLON_DATA_MORTON);
  {
    int d1 = x_axis == 0 ? x_length : (y_axis == 0 ? y_length : 1);
    int d2 = x_axis == 1 ? x_length : (y_axis == 1 ? y_length : 1);
    int d3 = x_axis == 2 ? x_length : (y_axis == 2 ? y_length : 1);
    /*  Fill with values 0, 1, 2 , ...  */
    t8_locidx_t         index;
    int x, y, z;
    for (x = 0; x < d1; ++x) {
      for (y = 0; y < d2; ++y) {
        for (z = 0; z < d3; ++z) {
          chunk->in[x][y][z][0] = (x * y_length + y) * dimension;
        }
      }
    }

    /* Print input data */
    t8_debugf ("Data:\n");
    for (x = 0; x < d1; ++x) {
      for (y = 0; y < d2; ++y) {
        for (z = 0; z < d3; ++z) {
          chunk->in[x][y][z][0] = (x * y_length + y) * dimension;
          t8_debugf ("[%d][%d][%d]: %.2f\n", x, y, z, chunk->in[x][y][z][0]);
        }
      }
    }

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

/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2024 the developers

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

#include "t8_vtk/t8_vtk_write_ASCII.hxx"
#include "t8_vtk/t8_vtk_writer_helper.hxx"
#include <t8_vtk.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_types/t8_vec.hxx>
#include "t8_forest/t8_forest_types.h"
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_cmesh/t8_cmesh_types.h"
#include <t8_schemes/t8_scheme.hxx>
#include <t8_geometry/t8_geometry_handler.hxx>

/* TODO: Currently we only use ASCII mode and no data compression.
 *       We also do not use sc_io to buffer our output stream. */

/* There are different cell data to write, e.g. connectivity, type, vertices, ...
 * The structure is always the same:
 * Iterate over the trees,
 *      iterate over the elements of that tree
 *          execute an element dependent part to write in the file.
 * In order to simplify writing this code, we put all the parts that are
 * repetitive in the function
 *  t8_forest_vtk_write_cell_data.
 * This function accepts a callback function, which is then executed for
 * each element. The callback function is defined below.
 */
/* TODO: As soon as we have element iterators we should restructure this concept
 * appropriately. */
typedef enum { T8_VTK_KERNEL_INIT, T8_VTK_KERNEL_EXECUTE, T8_VTK_KERNEL_CLEANUP } T8_VTK_KERNEL_MODUS;

/** Callback function prototype for writing cell data.
 * The function is executed for each element.
 * The callback can run in three different modi:
 *  INIT    - Called once, to (possibly) initialize the data pointer
 *  EXECUTE - Called for each element, the actual writing happens here.
 *  CLEANUP - Called once after all elements. Used to cleanup any memory
 *            allocated during INIT.
 * \param [in] forest The forest.
 * \param [in] ltree_id   A local treeid.
 * \param [in] tree   The local tree of the forest with id \a ltree_id.
 * \param [in] element_index An index of an element inside \a tree.
 * \param [in] element  A pointer to the current element.
 * \param [in] scheme       The eclass scheme of the current element.
 * \param [in] is_ghost Non-zero if the current element is a ghost element.
 *                      In this cas \a tree is NULL.
 *                      All ghost element will be traversed after all elements are
 * \param [in,out] vtufile The open file stream to which we write the forest.
 * \param [in,out] columns An integer counting the number of written columns.
 *                         The callback should increase this value by the number
 *                         of values written to the file.
 * \param [in,out] data    A pointer that the callback can modify at will.
 *                         Between modi INIT and CLEANUP, \a data will not be
 *                         modified outside of this callback.
 * \param [in]     modus   The modus in which the callback is called. See above.
 * \return                 True if successful, false if not (i.e. file i/o error).
 */
typedef int (*t8_forest_vtk_cell_data_kernel) (t8_forest_t forest, const t8_locidx_t ltree_id, const t8_tree_t tree,
                                               const t8_locidx_t element_index, const t8_element_t *element,
                                               const t8_eclass_t tree_class, const int is_ghost, FILE *vtufile,
                                               int *columns, void **data, T8_VTK_KERNEL_MODUS modus);

static t8_locidx_t
t8_forest_num_points (t8_forest_t forest, const int count_ghosts)
{
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  t8_locidx_t num_points = 0;

  for (t8_locidx_t itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    /* Get the tree that stores the elements */
    t8_tree_t tree = (t8_tree_t) t8_sc_array_index_locidx (forest->trees, itree);
    /* Get the scheme of the current tree */
    const size_t num_elements = t8_element_array_get_count (&tree->elements);
    for (t8_locidx_t ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
      const t8_element_t *elem = t8_element_array_index_locidx (&tree->elements, ielem);
      num_points += scheme->element_get_num_corners (tree_class, elem);
    }
  }
  if (count_ghosts) {
    T8_ASSERT (forest->ghosts != NULL);
    /* We also count the points of the ghost cells */
    const t8_locidx_t num_ghosts = t8_forest_ghost_num_trees (forest);
    for (t8_locidx_t itree = 0; itree < num_ghosts; itree++) {
      /* Get the element class of the ghost */
      const t8_eclass_t ghost_class = t8_forest_ghost_get_tree_class (forest, itree);
      const t8_element_array_t *ghost_elem = t8_forest_ghost_get_tree_elements (forest, itree);
      const size_t num_elements = t8_forest_ghost_tree_num_elements (forest, itree);
      for (t8_locidx_t ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
        const t8_element_t *elem = t8_element_array_index_locidx (ghost_elem, ielem);
        num_points += scheme->element_get_num_corners (ghost_class, elem);
      }
    }
  }
  return num_points;
}

static int
t8_forest_vtk_cells_vertices_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                     [[maybe_unused]] const t8_tree_t tree,
                                     [[maybe_unused]] const t8_locidx_t element_index, const t8_element_t *element,
                                     const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost, FILE *vtufile,
                                     int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  double element_coordinates[3];
  int num_el_vertices, ivertex;
  int freturn;
  t8_element_shape_t element_shape;

  if (modus != T8_VTK_KERNEL_EXECUTE) {
    /* Nothing to do if we are in Init or clean up mode */
    return 1;
  }

  /* TODO: be careful with pyramid class here.
   *       does this work too over tree->class or do we need something else?
   */
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  element_shape = scheme->element_get_shape (tree_class, element);
  num_el_vertices = t8_eclass_num_vertices[element_shape];
  for (ivertex = 0; ivertex < num_el_vertices; ivertex++) {
    const double *ref_coords = t8_forest_vtk_point_to_element_ref_coords[element_shape][ivertex];
    t8_forest_element_from_ref_coords (forest, ltree_id, element, ref_coords, 1, element_coordinates);
    freturn = fprintf (vtufile, "         ");
    if (freturn <= 0) {
      return 0;
    }
#ifdef T8_VTK_DOUBLES
    freturn = fprintf (vtufile, " %24.16e %24.16e %24.16e\n", element_coordinates[0], element_coordinates[1],
                       element_coordinates[2]);
#else
    freturn = fprintf (vtufile, " %16.8e %16.8e %16.8e\n", element_coordinates[0], element_coordinates[1],
                       element_coordinates[2]);
#endif
    if (freturn <= 0) {
      return 0;
    }
    /* We switch of the column control of the surrounding function
     * by keeping the columns value constant. */
    *columns = 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_connectivity_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                         [[maybe_unused]] const t8_tree_t tree,
                                         [[maybe_unused]] const t8_locidx_t element_index, const t8_element_t *element,
                                         const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost,
                                         FILE *vtufile, int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  int ivertex, num_vertices;
  int freturn;
  t8_locidx_t *count_vertices;
  t8_element_shape_t element_shape;

  if (modus == T8_VTK_KERNEL_INIT) {
    /* We use data to count the number of written vertices */
    *data = T8_ALLOC_ZERO (t8_locidx_t, 1);
    return 1;
  }
  else if (modus == T8_VTK_KERNEL_CLEANUP) {
    T8_FREE (*data);
    return 1;
  }
  T8_ASSERT (modus == T8_VTK_KERNEL_EXECUTE);

  count_vertices = (t8_locidx_t *) *data;
  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  element_shape = scheme->element_get_shape (tree_class, element);
  num_vertices = t8_eclass_num_vertices[element_shape];
  for (ivertex = 0; ivertex < num_vertices; ++ivertex, (*count_vertices)++) {
    freturn = fprintf (vtufile, " %ld", (long) *count_vertices);
    if (freturn <= 0) {
      return 0;
    }
  }
  *columns += t8_eclass_num_vertices[element_shape];
  return 1;
}

static int
t8_forest_vtk_cells_offset_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                   [[maybe_unused]] const t8_tree_t tree,
                                   [[maybe_unused]] const t8_locidx_t element_index, const t8_element_t *element,
                                   const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost, FILE *vtufile,
                                   int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  long long *offset;
  int freturn;
  int num_vertices;

  if (modus == T8_VTK_KERNEL_INIT) {
    *data = T8_ALLOC_ZERO (long long, 1);
    return true;
  }
  else if (modus == T8_VTK_KERNEL_CLEANUP) {
    T8_FREE (*data);
    return true;
  }
  T8_ASSERT (modus == T8_VTK_KERNEL_EXECUTE);

  offset = (long long *) *data;

  const t8_scheme *scheme = t8_forest_get_scheme (forest);
  num_vertices = t8_eclass_num_vertices[scheme->element_get_shape (tree_class, element)];
  *offset += num_vertices;
  freturn = fprintf (vtufile, " %lld", *offset);
  if (freturn <= 0) {
    return false;
  }
  *columns += 1;

  return 1;
}

static int
t8_forest_vtk_cells_type_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                 [[maybe_unused]] const t8_tree_t tree,
                                 [[maybe_unused]] const t8_locidx_t element_index, const t8_element_t *element,
                                 const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost, FILE *vtufile,
                                 int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  int freturn;
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    /* print the vtk type of the element */
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    freturn = fprintf (vtufile, " %d", t8_eclass_vtk_type[scheme->element_get_shape (tree_class, element)]);
    if (freturn <= 0) {
      return 0;
    }
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_level_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                  [[maybe_unused]] const t8_tree_t tree,
                                  [[maybe_unused]] const t8_locidx_t element_index, const t8_element_t *element,
                                  const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost, FILE *vtufile,
                                  int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    fprintf (vtufile, "%i ", scheme->element_get_level (tree_class, element));
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_rank_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                 [[maybe_unused]] const t8_tree_t tree,
                                 [[maybe_unused]] const t8_locidx_t element_index,
                                 [[maybe_unused]] const t8_element_t *element,
                                 [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] const int is_ghost,
                                 FILE *vtufile, int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    fprintf (vtufile, "%i ", forest->mpirank);
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_treeid_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                   [[maybe_unused]] const t8_tree_t tree,
                                   [[maybe_unused]] const t8_locidx_t element_index,
                                   [[maybe_unused]] const t8_element_t *element,
                                   [[maybe_unused]] const t8_eclass_t tree_class, const int is_ghost, FILE *vtufile,
                                   int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    long long tree_id;
    if (is_ghost) {
      /* For ghost elements we write -1 as the tree is */
      tree_id = -1;
    }
    else {
      /* Otherwise the global tree id */
      tree_id = (long long) ltree_id + forest->first_local_tree;
    }
    fprintf (vtufile, "%lli ", tree_id);
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_elementid_kernel (t8_forest_t forest, [[maybe_unused]] const t8_locidx_t ltree_id,
                                      const t8_tree_t tree, [[maybe_unused]] const t8_locidx_t element_index,
                                      [[maybe_unused]] const t8_element_t *element,
                                      [[maybe_unused]] const t8_eclass_t tree_class, const int is_ghost, FILE *vtufile,
                                      int *columns, [[maybe_unused]] void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    if (!is_ghost) {
      fprintf (vtufile, "%lli ",
               element_index + tree->elements_offset + (long long) t8_forest_get_first_local_element_id (forest));
    }
    else {
      fprintf (vtufile, "%lli ", (long long) -1);
    }
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_scalar_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                   [[maybe_unused]] const t8_tree_t tree, const t8_locidx_t element_index,
                                   [[maybe_unused]] const t8_element_t *element,
                                   [[maybe_unused]] const t8_eclass_t tree_class, const int is_ghost, FILE *vtufile,
                                   int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  double element_value = 0;
  t8_locidx_t scalar_index;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    /* For local elements access the data array, for ghosts, write 0 */
    if (!is_ghost) {
      scalar_index = t8_forest_get_tree_element_offset (forest, ltree_id) + element_index;
      element_value = ((double *) *data)[scalar_index];
    }
    else {
      element_value = 0;
    }
    fprintf (vtufile, "%g ", element_value);
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_vector_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                   [[maybe_unused]] const t8_tree_t tree, const t8_locidx_t element_index,
                                   [[maybe_unused]] const t8_element_t *element,
                                   [[maybe_unused]] const t8_eclass_t tree_class, const int is_ghost, FILE *vtufile,
                                   int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  double *element_values, null_vec[3] = { 0, 0, 0 };
  int dim, idim;
  t8_locidx_t tree_offset;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    dim = 3;
    T8_ASSERT (forest->dimension <= 3);
    /* For local elements access the data array, for ghosts, write 0 */
    if (!is_ghost) {
      tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);
      /* Get a pointer to the start of the element's vector data */
      element_values = ((double *) *data) + (tree_offset + element_index) * dim;
    }
    else {
      element_values = null_vec;
    }
    for (idim = 0; idim < dim; idim++) {
      fprintf (vtufile, "%g ", element_values[idim]);
    }
    *columns += dim;
  }
  return 1;
}

/* The point data version of the scalar kernel */
static int
t8_forest_vtk_vertices_scalar_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                      [[maybe_unused]] const t8_tree_t tree, const t8_locidx_t element_index,
                                      const t8_element_t *element, const t8_eclass_t tree_class, const int is_ghost,
                                      FILE *vtufile, int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  double element_value = 0;
  int num_vertex, ivertex;
  t8_locidx_t scalar_index;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    num_vertex = scheme->element_get_num_corners (tree_class, element);

    for (ivertex = 0; ivertex < num_vertex; ivertex++) {
      /* For local elements access the data array, for ghosts, write 0 */
      if (!is_ghost) {
        scalar_index = t8_forest_get_tree_element_offset (forest, ltree_id) + element_index;
        element_value = ((double *) *data)[scalar_index];
      }
      else {
        element_value = 0;
      }
      fprintf (vtufile, "%g ", element_value);
      *columns += 1;
    }
  }
  return 1;
}

/* The point data version of the vector kernel */
static int
t8_forest_vtk_vertices_vector_kernel (t8_forest_t forest, const t8_locidx_t ltree_id,
                                      [[maybe_unused]] const t8_tree_t tree, const t8_locidx_t element_index,
                                      const t8_element_t *element, const t8_eclass_t tree_class, const int is_ghost,
                                      FILE *vtufile, int *columns, void **data, T8_VTK_KERNEL_MODUS modus)
{
  double *element_values, null_vec[3] = { 0, 0, 0 };
  int dim, idim;
  int num_vertex, ivertex;
  t8_locidx_t tree_offset;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    const t8_scheme *scheme = t8_forest_get_scheme (forest);
    num_vertex = scheme->element_get_num_corners (tree_class, element);
    for (ivertex = 0; ivertex < num_vertex; ivertex++) {
      dim = 3;
      T8_ASSERT (forest->dimension <= 3);
      /* For local elements access the data array, for ghosts, write 0 */
      if (!is_ghost) {
        tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);
        /* Get a pointer to the start of the element's vector data */
        element_values = ((double *) *data) + (tree_offset + element_index) * dim;
      }
      else {
        element_values = null_vec;
      }
      for (idim = 0; idim < dim; idim++) {
        fprintf (vtufile, "%g ", element_values[idim]);
      }
      *columns += dim;
    }
  }
  return 1;
}

/* Iterate over all cells and write cell data to the file using
 * the cell_data_kernel as callback */
static int
t8_forest_vtk_write_cell_data (t8_forest_t forest, FILE *vtufile, const char *dataname, const char *datatype,
                               const char *component_string, const int max_columns,
                               t8_forest_vtk_cell_data_kernel kernel, const int write_ghosts, void *udata)
{
  int freturn;
  int countcols;
  t8_tree_t tree;
  t8_locidx_t itree, ighost;
  t8_locidx_t element_index, elems_in_tree;
  t8_locidx_t num_local_trees, num_ghost_trees;
  t8_element_t *element;
  void *data = NULL;

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = fprintf (vtufile,
                     "        <DataArray type=\"%s\" "
                     "Name=\"%s\" %s format=\"ascii\">\n         ",
                     datatype, dataname, component_string);
  if (freturn <= 0) {
    return 0;
  }

  /* if udata != NULL, use it as the data pointer, in this case, the kernel
   * should not modify it */
  if (udata != NULL) {
    data = udata;
  }

  /* Call the kernel in initialization modus to possibly initialize the
   * data pointer */
  kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_COUNT, 0, NULL, NULL, &data, T8_VTK_KERNEL_INIT);
  /* We iterate over the trees and count each trees vertices,
   * we add this to the already counted vertices and write it to the file */
  /* TODO: replace with an element iterator */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, countcols = 0; itree < num_local_trees; itree++) {
    /* Get the tree that stores the elements */
    tree = t8_forest_get_tree (forest, itree);
    /* Get the eclass scheme of the tree */
    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
    elems_in_tree = (t8_locidx_t) t8_element_array_get_count (&tree->elements);
    for (element_index = 0; element_index < elems_in_tree; element_index++) {
      /* Get a pointer to the element */
      element = t8_forest_get_element (forest, tree->elements_offset + element_index, NULL);
      T8_ASSERT (element != NULL);
      /* Execute the given callback on each element */
      if (!kernel (forest, itree, tree, element_index, element, tree_class, 0, vtufile, &countcols, &data,
                   T8_VTK_KERNEL_EXECUTE)) {
        /* call the kernel in clean-up modus */
        kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_COUNT, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
        return 0;
      }
      /* After max_columns we break the line */
      if (!(countcols % max_columns)) {
        freturn = fprintf (vtufile, "\n         ");
        if (freturn <= 0) {
          /* call the kernel in clean-up modus */
          kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_COUNT, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
          return 0;
        }
      }
    } /* element loop ends here */
    if (freturn <= 0) {
      /* call the kernel in clean-up modus */
      kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_INVALID, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
      return 0;
    }
  } /* tree loop ends here */

  if (write_ghosts) {
    t8_locidx_t num_ghosts_in_tree;
    /* Iterate over the ghost elements */
    /* TODO: replace with an element iterator */
    num_ghost_trees = t8_forest_ghost_num_trees (forest);
    for (ighost = 0; ighost < num_ghost_trees; ighost++) {
      /* Get the eclass of the ghost tree */
      const t8_eclass_t ghost_eclass = t8_forest_ghost_get_tree_class (forest, ighost);
      /* The number of ghosts in this tree */
      num_ghosts_in_tree = t8_forest_ghost_tree_num_elements (forest, ighost);
      for (element_index = 0; element_index < num_ghosts_in_tree; element_index++) {
        /* Get a pointer to the element */
        element = t8_forest_ghost_get_element (forest, ighost, element_index);
        /* Execute the given callback on each element */
        if (!kernel (forest, ighost + num_local_trees, NULL, element_index, element, ghost_eclass, 1, vtufile,
                     &countcols, &data, T8_VTK_KERNEL_EXECUTE)) {
          /* call the kernel in clean-up modus */
          kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_INVALID, 1, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
          return 0;
        }
        /* After max_columns we break the line */
        if (!(countcols % max_columns)) {
          freturn = fprintf (vtufile, "\n         ");
          if (freturn <= 0) {
            /* call the kernel in clean-up modus */
            kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_INVALID, 1, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
            return 0;
          }
        }
      } /* element loop ends here */
      if (freturn <= 0) {
        /* call the kernel in clean-up modus */
        kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_INVALID, 1, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
        return 0;
      }
    } /* ghost loop ends here */
  }   /* write_ghosts ends here */
  /* call the kernel in clean-up modus */
  kernel (NULL, 0, NULL, 0, NULL, T8_ECLASS_INVALID, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
  freturn = fprintf (vtufile, "\n        </DataArray>\n");
  if (freturn <= 0) {
    return 0;
  }

  return 1;
}

/* Write the cell data to an open file stream.
 * Returns true on success and zero otherwise.
 * After completion the file will remain open, whether writing
 * cells was successful or not. */
static int
t8_forest_vtk_write_cells (t8_forest_t forest, FILE *vtufile, const int write_treeid, const int write_mpirank,
                           const int write_level, const int write_element_id, const int write_ghosts,
                           const int num_data, t8_vtk_data_field_t *data)
{
  int freturn;
  int idata;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (vtufile != NULL);

  freturn = fprintf (vtufile, "      <Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "connectivity", T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_connectivity_kernel, write_ghosts, NULL);
  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the connectivity */

  /* Write the offsets, that is for each tree the index of the first entry
   * in the connectivity output that
   * does not refer to a vertex of the tree anymore.
   * For example if the trees are a square and a triangle, the offsets would
   * be 4 and 7, since indices 0,1,2,3 refer to the vertices of the square
   * and indices 4,5,6 to the indices of the triangle. */
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "offsets", T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_offset_kernel, write_ghosts, NULL);
  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the offsets */

  /* Write the element types. The type specifies the element class, thus
   * square/triangle/tet etc. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "types", "Int32", "", 8, t8_forest_vtk_cells_type_kernel,
                                           write_ghosts, NULL);

  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the types */
  freturn = fprintf (vtufile, "      </Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* clang-format off */
  freturn = fprintf (vtufile, "      <CellData Scalars =\"%s%s\">\n", "treeid,mpirank,level",
                     (write_element_id ? "id" : ""));
  /* clang-format on */
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  if (write_treeid) {
    /* Write the tree ids. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "treeid", T8_VTK_GLOIDX, "", 8,
                                             t8_forest_vtk_cells_treeid_kernel, write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the tree ids */
  }
  if (write_mpirank) {
    /* Write the mpiranks. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "mpirank", "Int32", "", 8,
                                             t8_forest_vtk_cells_rank_kernel, write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the mpiranks */
  }
  if (write_level) {
    /* Write the element refinement levels. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "level", "Int32", "", 8, t8_forest_vtk_cells_level_kernel,
                                             write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the levels */
  }

  if (write_element_id) {
    /* Write the element ids. */
    const char *datatype;

    /* Use 32 bit ints if the global element count fits, 64 bit otherwise. */
    datatype = forest->global_num_elements > T8_LOCIDX_MAX ? T8_VTK_GLOIDX : T8_VTK_LOCIDX;
    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "element_id", datatype, "", 8,
                                             t8_forest_vtk_cells_elementid_kernel, write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }

    /* Done with writing the element ids */
  }
  /* Write the user defined data fields per element */
  for (idata = 0; idata < num_data; idata++) {
    if (data[idata].type == T8_VTK_SCALAR) {
      freturn = t8_forest_vtk_write_cell_data (forest, vtufile, data[idata].description, T8_VTK_FLOAT_NAME, "", 8,
                                               t8_forest_vtk_cells_scalar_kernel, write_ghosts, data[idata].data);
    }
    else {
      char component_string[BUFSIZ];
      T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
      snprintf (component_string, BUFSIZ, "NumberOfComponents=\"3\"");
      freturn = t8_forest_vtk_write_cell_data (forest, vtufile, data[idata].description, T8_VTK_FLOAT_NAME,
                                               component_string, 8 * forest->dimension,
                                               t8_forest_vtk_cells_vector_kernel, write_ghosts, data[idata].data);
    }
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
  }

  freturn = fprintf (vtufile, "      </CellData>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  /* Function completed successfully */
  return 1;
t8_forest_vtk_cell_failure:
  /* Something went wrong */
  t8_errorf ("Error when writing cell data to forest vtk file.\n");
  return 0;
}

/* Write the cell data to an open file stream.
 * Returns true on success and zero otherwise.
 * After completion the file will remain open, whether writing
 * cells was successful or not. */
static int
t8_forest_vtk_write_points (t8_forest_t forest, FILE *vtufile, const int write_ghosts, const int num_data,
                            t8_vtk_data_field_t *data)
{
  int freturn;
  int sreturn;
  int idata;
  char description[BUFSIZ];

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (vtufile != NULL);

  /* Write the vertex coordinates */

  freturn = fprintf (vtufile, "      <Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "Position", T8_VTK_FLOAT_NAME, "NumberOfComponents=\"3\"",
                                           8, t8_forest_vtk_cells_vertices_kernel, write_ghosts, NULL);
  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  freturn = fprintf (vtufile, "      </Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done writing vertex coordinates */

  /* Write the user defined data fields per element */
  if (num_data > 0) {
    freturn = fprintf (vtufile, "      <PointData>\n");
    for (idata = 0; idata < num_data; idata++) {
      if (data[idata].type == T8_VTK_SCALAR) {
        /* Write the description string. */
        sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
        }
        freturn = t8_forest_vtk_write_cell_data (forest, vtufile, description, T8_VTK_FLOAT_NAME, "", 8,
                                                 t8_forest_vtk_vertices_scalar_kernel, write_ghosts, data[idata].data);
      }
      else {
        char component_string[BUFSIZ];
        T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
        snprintf (component_string, BUFSIZ, "NumberOfComponents=\"3\"");
        /* Write the description string. */
        sreturn = snprintf (description, BUFSIZ, "%s_%s", data[idata].description, "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          /* Note: gcc >= 7.1 prints a warning if we 
           * do not check the return value of snprintf. */
          t8_debugf ("Warning: Truncated vtk point data description to '%s'\n", description);
        }

        freturn = t8_forest_vtk_write_cell_data (forest, vtufile, description, T8_VTK_FLOAT_NAME, component_string,
                                                 8 * forest->dimension, t8_forest_vtk_vertices_vector_kernel,
                                                 write_ghosts, data[idata].data);
      }
      if (!freturn) {
        goto t8_forest_vtk_cell_failure;
      }
    }
    freturn = fprintf (vtufile, "      </PointData>\n");
  }
  /* Function completed successfully */
  return 1;
t8_forest_vtk_cell_failure:
  /* Something went wrong */
  t8_errorf ("Error when writing cell data to forest vtk file.\n");
  return 0;
}

int
t8_forest_vtk_write_ASCII (t8_forest_t forest, const char *fileprefix, const int write_treeid, const int write_mpirank,
                           const int write_level, const int write_element_id, int write_ghosts, const int num_data,
                           t8_vtk_data_field_t *data)
{
  FILE *vtufile = NULL;
  t8_locidx_t num_elements, num_points;
  char vtufilename[BUFSIZ];
  int freturn;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (fileprefix != NULL);
  if (forest->cmesh->geometry_handler == nullptr) {
    t8_errorf ("Error: Forest could not be written as vtk since it does not have geometries.\n");
    return false;
  }
  if (forest->cmesh->geometry_handler->get_num_geometries () == 0) {
    t8_errorf ("Error: Forest could not be written as vtk since it does not have geometries.\n");
    return false;
  }
  if (forest->ghosts == NULL || forest->ghosts->num_ghosts_elements == 0) {
    /* Never write ghost elements if there aren't any */
    write_ghosts = 0;
  }
  T8_ASSERT (forest->ghosts != NULL || !write_ghosts);

  /* process 0 creates the .pvtu file */
  if (forest->mpirank == 0) {
    if (t8_write_pvtu (fileprefix, forest->mpisize, write_treeid, write_mpirank, write_level, write_element_id,
                       num_data, data)) {
      t8_errorf ("Error when writing file %s.pvtu\n", fileprefix);
      goto t8_forest_vtk_failure;
    }
  }

  /* The local number of elements */
  num_elements = t8_forest_get_local_num_elements (forest);
  if (write_ghosts) {
    num_elements += t8_forest_get_num_ghosts (forest);
  }
  /* The local number of points, counted with multiplicity */
  num_points = t8_forest_num_points (forest, write_ghosts);

  /* The filename for this processes file */
  freturn = snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", fileprefix, forest->mpirank);
  if (freturn >= BUFSIZ) {
    t8_errorf ("Error when writing vtu file. Filename too long.\n");
    goto t8_forest_vtk_failure;
  }

  /* Open the vtufile to write to */
  vtufile = fopen (vtufilename, "w");
  if (vtufile == NULL) {
    t8_errorf ("Error when opening file %s\n", vtufilename);
    goto t8_forest_vtk_failure;
  }
  /* Write the header information in the .vtu file.
   * xml type, Unstructured grid and number of points and elements. */
  freturn = fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }
  freturn = fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }
#ifdef SC_IS_BIGENDIAN
  freturn = fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
  freturn = fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }
  freturn = fprintf (vtufile, "  <UnstructuredGrid>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }
  freturn = fprintf (vtufile, "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) num_points,
                     (long long) num_elements);
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }
  /* write the point data */
  if (!t8_forest_vtk_write_points (forest, vtufile, write_ghosts, num_data, data)) {
    /* writings points was not successful */
    goto t8_forest_vtk_failure;
  }
  /* write the cell data */
  if (!t8_forest_vtk_write_cells (forest, vtufile, write_treeid, write_mpirank, write_level, write_element_id,
                                  write_ghosts, num_data, data)) {
    /* Writing cells was not successful */
    goto t8_forest_vtk_failure;
  }

  freturn = fprintf (vtufile, "    </Piece>\n"
                              "  </UnstructuredGrid>\n"
                              "</VTKFile>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }

  freturn = fclose (vtufile);
  /* We set it not NULL, even if fclose was not successful, since then any
   * following call to fclose would result in undefined behaviour. */
  vtufile = NULL;
  if (freturn != 0) {
    /* Closing failed, this usually means that the final write operation could
     * not be completed. */
    t8_global_errorf ("Error when closing file %s\n", vtufilename);
    goto t8_forest_vtk_failure;
  }
  /* Writing was successful */
  return 1;
t8_forest_vtk_failure:
  if (vtufile != NULL) {
    fclose (vtufile);
  }
  t8_errorf ("Error when writing vtk file.\n");
  return 0;
}

/* Return the local number of vertices in a cmesh.
 * \param [in] cmesh       The cmesh to be considered.
 * \param [in] count_ghosts If true, we also count the vertices of the ghost trees.
 * \return                 The number of vertices associated to \a cmesh.
 * \a cmesh must be committed before calling this function.
 */
static t8_gloidx_t
t8_cmesh_get_num_vertices (const t8_cmesh_t cmesh, const int count_ghosts)
{
  int iclass;
  t8_eclass_t ghost_class;
  t8_gloidx_t num_vertices = 0;
  t8_locidx_t ighost;
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (cmesh->committed);

  for (iclass = T8_ECLASS_ZERO; iclass < T8_ECLASS_COUNT; iclass++) {
    num_vertices += t8_eclass_num_vertices[iclass] * cmesh->num_local_trees_per_eclass[iclass];
  }
  if (count_ghosts) {
    /* Also count the vertices of the ghost trees */
    for (ighost = 0; ighost < t8_cmesh_get_num_ghosts (cmesh); ighost++) {
      ghost_class = t8_cmesh_get_ghost_class (cmesh, ighost);
      num_vertices += t8_eclass_num_vertices[ghost_class];
    }
  }

  return num_vertices;
}

static int
t8_cmesh_vtk_write_file_ext (const t8_cmesh_t cmesh, const char *fileprefix, const int write_ghosts)
{
  T8_ASSERT (cmesh != NULL);
  T8_ASSERT (t8_cmesh_is_committed (cmesh));
  T8_ASSERT (fileprefix != NULL);

  if (cmesh->geometry_handler == nullptr) {
    t8_errorf ("Error: Cmesh could not be written as vtk since it does not have geometries.\n");
    return false;
  }
  if (cmesh->geometry_handler->get_num_geometries () == 0) {
    t8_errorf ("Error: Cmesh could not be written as vtk since it does not have geometries.\n");
    return false;
  }

  if (cmesh->mpirank == 0) {
    /* Write the pvtu header file. */
    int num_ranks_that_write = cmesh->set_partition ? cmesh->mpisize : 1;
    if (t8_write_pvtu (fileprefix, num_ranks_that_write, 1, 1, 0, 0, 0, NULL)) {
      SC_ABORTF ("Error when writing file %s.pvtu\n", fileprefix);
    }
  }
  /* If the cmesh is replicated only rank 0 prints it,
   * otherwise each process prints its part of the cmesh.*/
  if (cmesh->mpirank == 0 || cmesh->set_partition) {
    char vtufilename[BUFSIZ];
    FILE *vtufile;
    t8_locidx_t num_vertices, ivertex;
    t8_locidx_t num_trees;
    t8_ctree_t tree;
    double x, y, z;
    double *vertices, *vertex;
    int k, sk;
    long long offset, count_vertices;
    t8_locidx_t ighost, num_ghosts = 0, num_loc_trees;
#ifdef T8_ENABLE_DEBUG
    t8_cghost_t ghost;
#endif
    t8_eclass_t eclass;

    num_vertices = t8_cmesh_get_num_vertices (cmesh, write_ghosts);
    num_trees = t8_cmesh_get_num_local_trees (cmesh);
    if (write_ghosts) {
      num_trees += t8_cmesh_get_num_ghosts (cmesh);
    }

    snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", fileprefix, cmesh->mpirank);
    vtufile = fopen (vtufilename, "wb");
    if (vtufile == NULL) {
      t8_global_errorf ("Could not open file %s for output.\n", vtufilename);
      return 0;
    }
    fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#ifdef SC_IS_BIGENDIAN
    fprintf (vtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
#endif
    fprintf (vtufile, "  <UnstructuredGrid>\n");
    fprintf (vtufile, "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) num_vertices,
             (long long) num_trees);
    fprintf (vtufile, "      <Points>\n");

    /* write point position data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\">\n",
             T8_VTK_FLOAT_NAME, T8_VTK_FORMAT_STRING);

    for (tree = t8_cmesh_get_first_tree (cmesh); tree != NULL; tree = t8_cmesh_get_next_tree (cmesh, tree)) {
      /*  TODO: Use new geometry here. Need cmesh_get_reference coords function. */
      vertices = t8_cmesh_get_tree_vertices (cmesh, tree->treeid);
      for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass]; ivertex++) {
        vertex = vertices + 3 * t8_eclass_t8_to_vtk_corner_number[tree->eclass][ivertex];
        x = vertex[0];
        y = vertex[1];
        z = vertex[2];
#ifdef T8_VTK_DOUBLES
        fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
        fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
      }
    } /* end tree loop */
    if (write_ghosts) {

      /* Write the vertices of the ghost trees */
      num_ghosts = t8_cmesh_get_num_ghosts (cmesh);
      num_loc_trees = t8_cmesh_get_num_local_trees (cmesh);
      for (ighost = 0; ighost < num_ghosts; ighost++) {
        /* Get the eclass of this ghost */
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        /* Get a pointer to this ghosts vertices */
        vertices = (double *) t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0, ighost + num_loc_trees);
        T8_ASSERT (vertices != NULL);
        /* TODO: This code is duplicated above */
        for (ivertex = 0; ivertex < t8_eclass_num_vertices[eclass]; ivertex++) {
          vertex = vertices + 3 * t8_eclass_vtk_to_t8_corner_number[eclass][ivertex];
          x = vertex[0];
          y = vertex[1];
          z = vertex[2];
#ifdef T8_VTK_DOUBLES
          fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
          fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
        }
      } /* end ghost loop */
    }
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Points>\n");
    fprintf (vtufile, "      <Cells>\n");

    /* write connectivity data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"connectivity\""
             " format=\"%s\">\n",
             T8_VTK_LOCIDX, T8_VTK_FORMAT_STRING);
    for (tree = t8_cmesh_get_first_tree (cmesh), count_vertices = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree)) {
      fprintf (vtufile, "         ");
      for (k = 0; k < t8_eclass_num_vertices[tree->eclass]; ++k, count_vertices++) {
        fprintf (vtufile, " %lld", count_vertices);
      }
      fprintf (vtufile, "\n");
    }
    if (write_ghosts) {
      /* Write the ghost connectivity */
      for (ighost = 0; ighost < num_ghosts; ighost++) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        fprintf (vtufile, "         ");
        for (k = 0; k < t8_eclass_num_vertices[eclass]; ++k, count_vertices++) {
          fprintf (vtufile, " %lld", count_vertices);
        }
        fprintf (vtufile, "\n");
      }
    }
    fprintf (vtufile, "        </DataArray>\n");

    /* write offset data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"offsets\""
             " format=\"%s\">\n",
             T8_VTK_LOCIDX, T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      offset += t8_eclass_num_vertices[tree->eclass];
      fprintf (vtufile, " %lld", offset);
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset data */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        offset += t8_eclass_num_vertices[eclass];
        fprintf (vtufile, " %lld", offset);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    /* write type data */
    fprintf (vtufile,
             "        <DataArray type=\"UInt8\" Name=\"types\""
             " format=\"%s\">\n",
             T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      fprintf (vtufile, " %d", t8_eclass_vtk_type[tree->eclass]);
      if (!(sk % 20) && tree->treeid != (cmesh->num_local_trees - 1))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset types */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        eclass = t8_cmesh_get_ghost_class (cmesh, ighost);
        fprintf (vtufile, " %d", t8_eclass_vtk_type[eclass]);
        if (!(sk % 20) && ighost != (num_ghosts - 1))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </Cells>\n");
    /* write treeif data */
    fprintf (vtufile, "      <CellData Scalars=\"treeid,mpirank\">\n");
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n",
             T8_VTK_GLOIDX, T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      /* Since tree_id is actually 64 Bit but we store it as 32, we have to check
       * that we do not get into conversion errors */
      /* TODO: We switched to 32 Bit because Paraview could not handle 64 well enough.
       */
      T8_ASSERT (tree->treeid + cmesh->first_tree == (t8_gloidx_t) ((long) tree->treeid + cmesh->first_tree));
      fprintf (vtufile, " %ld", static_cast<long> (tree->treeid + cmesh->first_tree));
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* ghost offset types */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
#ifdef T8_ENABLE_DEBUG
        ghost = t8_cmesh_trees_get_ghost (cmesh->trees, ighost);
        /* Check for conversion errors */
        T8_ASSERT (ghost->treeid == (t8_gloidx_t) ((long) ghost->treeid));
#endif
        /* Write -1 as tree_id so that we can distinguish ghosts from normal trees
         * in the vtk file */
        fprintf (vtufile, " %ld", (long) -1);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    /* write mpirank data */
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n",
             "Int32", T8_VTK_FORMAT_STRING);
    fprintf (vtufile, "         ");
    for (tree = t8_cmesh_get_first_tree (cmesh), sk = 1, offset = 0; tree != NULL;
         tree = t8_cmesh_get_next_tree (cmesh, tree), ++sk) {
      fprintf (vtufile, " %i", cmesh->mpirank);
      if (!(sk % 8))
        fprintf (vtufile, "\n         ");
    }
    if (write_ghosts) {
      /* write our rank for each ghost */
      for (ighost = 0; ighost < num_ghosts; ighost++, ++sk) {
        fprintf (vtufile, " %i", cmesh->mpirank);
        if (!(sk % 8))
          fprintf (vtufile, "\n         ");
      }
    }
    fprintf (vtufile, "\n");
    fprintf (vtufile, "        </DataArray>\n");
    fprintf (vtufile, "      </CellData>\n");
    /* write type data */
    fprintf (vtufile, "    </Piece>\n");
    fprintf (vtufile, "  </UnstructuredGrid>\n");
    fprintf (vtufile, "</VTKFile>\n");
    fclose (vtufile);
  }
  return 1;
}

int
t8_cmesh_vtk_write_ASCII (const t8_cmesh_t cmesh, const char *fileprefix)
{
  return t8_cmesh_vtk_write_file_ext (cmesh, fileprefix, 1);
}

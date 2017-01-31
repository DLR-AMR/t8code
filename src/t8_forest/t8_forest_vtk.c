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

#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_forest_types.h"

/* TODO: Currently we only use ASCII mode and no data compression.
 *       We also do not use sc_io to buffer our output stream. */

/* There are different cell data to write, e.g. connectivity, type, vertices, ...
 * The structure is always the same:
 * Iterate over the trees,
 *      iterate over the elements of that tree
 *          execute an element dependent part to write in the file.
 * In order to simplify writing this code, we put all the parts that are
 * repetetive in the function
 *  t8_forest_vtk_write_cell_data.
 * This function accepts a callback function, which is then executed for
 * each element. The callback function is defined below.
 */
/* TODO: As soon as we have element iterators we should restructure this concept
 * appropiately. */
typedef enum
{
  T8_VTK_KERNEL_INIT,
  T8_VTK_KERNEL_EXECUTE,
  T8_VTK_KERNEL_CLEANUP
} T8_VTK_KERNEL_MODUS;

/* Callback function prototype for writing cell data.
 * The function is executed for each element.
 * The callback can run in three different modi:
 *  INIT    - Called once, to (possibly) initialize the data pointer
 *  EXECUTE - Called for each element, the actual writing happens here.
 *  CLEAUP  - Called once after all elements. Used to cleanup any memory
 *            allocated during INIT.
 * \param [in] forest The forest.
 * \param [in] ltree_id   A local treeid.
 * \param [in] tree   The local tree of the forest with id \a ltree_id.
 * \param [in] element_index An index of an element inside \a tree.
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
typedef int         (*t8_forest_vtk_cell_data_kernel) (t8_forest_t forest,
                                                       t8_locidx_t ltree_id,
                                                       t8_tree_t tree,
                                                       t8_locidx_t
                                                       element_index,
                                                       FILE * vtufile,
                                                       int *columns,
                                                       void **data,
                                                       T8_VTK_KERNEL_MODUS
                                                       modus);

static              t8_locidx_t
t8_forest_num_points (t8_forest_t forest)
{
  t8_locidx_t         itree, num_points;
  t8_tree_t           tree;

  num_points = 0;
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* Get the tree that stores the elements */
    tree = (t8_tree_t) t8_sc_array_index_topidx (forest->trees, itree);
    /* TODO: This will cause problems when pyramids are introduced. */
    num_points += t8_eclass_num_vertices[tree->eclass] *
      tree->elements.elem_count;
  }
  return num_points;
}

static int
t8_forest_vtk_cells_vertices_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                     t8_tree_t tree,
                                     t8_locidx_t element_index,
                                     FILE * vtufile, int *columns,
                                     void **data, T8_VTK_KERNEL_MODUS modus)
{
  struct t8_forest_vtk_vertices_t
  {
    t8_tree_t           tree;   /* Store the last tree with which the kernel was called */
    double              tree_vertices[T8_ECLASS_MAX_CORNERS * 3];       /* Stores the vertex coordinates of the tree */
  }                  *vertex_data;

  double              element_coordinates[3];
  int                 num_tree_vertices, ivertex;
  t8_element_t       *element;
  int                 freturn;

  if (modus == T8_VTK_KERNEL_INIT) {
    /* We initialize the user data to store NULL as the current tree */
    *data = T8_ALLOC_ZERO (struct t8_forest_vtk_vertices_t, 1);
    return 1;
  }
  else if (modus == T8_VTK_KERNEL_CLEANUP) {
    T8_FREE (*data);
    return 1;
  }

  T8_ASSERT (modus == T8_VTK_KERNEL_EXECUTE);
  vertex_data = (struct t8_forest_vtk_vertices_t *) *data;
  if (tree != vertex_data->tree) {
    t8_ctree_t          ctree;
    t8_cmesh_t          cmesh;
    double             *temp_vertices;
    /* The current tree is not the tree that we stored from
     * the last call to this function */
    vertex_data->tree = tree;
    /* get the coarse mesh tree */
    cmesh = t8_forest_get_cmesh (forest);
    ctree = t8_cmesh_get_tree (cmesh,
                               t8_forest_ltreeid_to_cmesh_ltreeid (forest,
                                                                   ltree_id));
    /* Get corner coordinates of tree */
    /* *INDENT-OFF* */
    /* indent bug */
    temp_vertices = ((double *)
                   t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0,
                                           ctree->treeid));
    /* *INDENT-ON* */
    /* Copy the tree's vertex coordinates into the struct of the data pointer */
    num_tree_vertices = t8_eclass_num_vertices[tree->eclass];
    memcpy (vertex_data->tree_vertices, temp_vertices, sizeof (*temp_vertices)
            * num_tree_vertices * 3);
  }

  element = (t8_element_t *) sc_array_index (&tree->elements, element_index);
  /* TODO: be careful with pyramid class here.
   *       does this work too over tree->class or do we need something else?
   */

  for (ivertex = 0;
       ivertex < t8_eclass_num_vertices[vertex_data->tree->eclass];
       ivertex++) {
    t8_forest_element_coordinate (forest, ltree_id, element,
                                  vertex_data->tree_vertices,
                                  t8_eclass_vtk_corner_number
                                  [vertex_data->tree->eclass]
                                  [ivertex], element_coordinates);
    freturn = fprintf (vtufile, "         ");
    if (freturn <= 0) {
      return 0;
    }
#ifdef T8_VTK_DOUBLES
    freturn = fprintf (vtufile, " %24.16e %24.16e %24.16e\n",
                       element_coordinates[0], element_coordinates[1],
                       element_coordinates[2]);
#else
    freturn = fprintf (vtufile, " %16.8e %16.8e %16.8e\n",
                       element_coordinates[0], element_coordinates[1],
                       element_coordinates[2]);
#endif
    if (freturn <= 0) {
      return 0;
    }
    /* We switch of the colum control of the surrounding function
     * by keeping the columns value constant. */
    *columns = 1;
  }
  return 1;
}

#if 0
/* Write vertex coordinates into the already opened file.
 * Returns true when successful, false otherwise.
 * After completion the file will remain open, whether writing
 * vertices was successful or not. */
static int
t8_forest_vtk_write_vertices (t8_forest_t forest, FILE * vtufile)
{
  t8_element_t       *ielement;
  t8_tree_t           tree;
  t8_locidx_t         itree, ivertex;
  t8_locidx_t         element_index;
  t8_cmesh_t          cmesh;
  t8_ctree_t          ctree;
  double             *vertices, coordinates[3];
  double              x, y, z;
  int                 freturn;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (vtufile != NULL);
  cmesh = forest->cmesh;
  freturn = fprintf (vtufile, "      <Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_vertex_failure;
  }

  /* write point position data */
  freturn =
    fprintf (vtufile,
             "        <DataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\">\n", T8_VTK_FLOAT_NAME,
             T8_VTK_FORMAT_STRING);
  if (freturn <= 0) {
    goto t8_forest_vtk_vertex_failure;
  }
  /* To get the point position data, we iterate over each tree and
   * over each element in this tree. For each element we compute
   * the coordinates of its corner vertices */
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* get the coarse mesh tree */
    ctree = t8_cmesh_get_tree (cmesh,
                               t8_forest_ltreeid_to_cmesh_ltreeid (forest,
                                                                   itree));
    /* Get corner coordinates of tree */
    /* *INDENT-OFF* */
    /* indent bug */
    vertices = ((double *)
                t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0,
                                        ctree->treeid));
    /* *INDENT-ON* */
    /* Get the tree that stores the elements */
    tree = t8_forest_get_tree (forest, itree);
    /* Check whether an element exist and then get the first one */
    /* TODO: use an element iterator here! */
    if (tree->elements.elem_count > 0) {
      ielement = (t8_element_t *) sc_array_index (&tree->elements, 0);
    }
    else {
      ielement = NULL;
    }
    element_index = 0;
    while (ielement != NULL) {
      /* TODO: be careful with pyramid class here.
       *       does this work too over tree->class or do we need something else?
       */
      for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass];
           ivertex++) {
        t8_forest_element_coordinate (forest, itree, ielement,
                                      vertices,
                                      t8_eclass_vtk_corner_number
                                      [tree->eclass]
                                      [ivertex], coordinates);
        x = coordinates[0];
        y = coordinates[1];
        z = coordinates[2];
#ifdef T8_VTK_DOUBLES
        freturn =
          fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", x, y, z);
#else
        freturn =
          fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", x, y, z);
#endif
        if (freturn <= 0) {
          goto t8_forest_vtk_vertex_failure;
        }
      }
      element_index++;
      ielement =
        element_index >=
        (t8_locidx_t) tree->elements.elem_count ? NULL : (t8_element_t *)
        t8_sc_array_index_locidx (&tree->elements, element_index);
    }
    /* loop over tree ends here */
  }
  freturn = fprintf (vtufile, "      </Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_vertex_failure;
  }

  /* Function completed successfuly */
  return 1;

t8_forest_vtk_vertex_failure:
  /* Something went wrong */
  t8_errorf ("Error when writing vertices to forest vtk file.\n");
  return 0;
}
#endif

static int
t8_forest_vtk_cells_connectivity_kernel (t8_forest_t forest,
                                         t8_locidx_t ltree_id, t8_tree_t tree,
                                         t8_locidx_t element_index,
                                         FILE * vtufile, int *columns,
                                         void **data,
                                         T8_VTK_KERNEL_MODUS modus)
{
  int                 ivertex;
  int                 freturn;
  t8_locidx_t        *count_vertices;

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

  count_vertices = (t8_locidx_t *) * data;

  /* TODO: This will definitely break with pyramids */
  for (ivertex = 0; ivertex < t8_eclass_num_vertices[tree->eclass];
       ++ivertex, (*count_vertices)++) {
    freturn = fprintf (vtufile, " %ld", (long) *count_vertices);
    if (freturn <= 0) {
      return 0;
    }
  }
  *columns += t8_eclass_num_vertices[tree->eclass];
  return 1;
}

static int
t8_forest_vtk_cells_offset_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                   t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  long long          *offset;
  int                 freturn;

  if (modus == T8_VTK_KERNEL_INIT) {
    *data = T8_ALLOC_ZERO (long long, 1);
    return 1;
  }
  else if (modus == T8_VTK_KERNEL_CLEANUP) {
    T8_FREE (*data);
    return 1;
  }
  T8_ASSERT (modus == T8_VTK_KERNEL_EXECUTE);

  offset = (long long *) *data;

  /* TODO: This will also break with pyramids! */
  *offset += t8_eclass_num_vertices[tree->eclass];
  freturn = fprintf (vtufile, " %lld", *offset);
  if (freturn <= 0) {
    return 0;
  }
  *columns += 1;

  return 1;
}

static int
t8_forest_vtk_cells_type_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                 t8_tree_t tree,
                                 t8_locidx_t element_index,
                                 FILE * vtufile, int *columns,
                                 void **data, T8_VTK_KERNEL_MODUS modus)
{
  int                 freturn;
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    /* print the vtk type of the element */
    freturn = fprintf (vtufile, " %d", t8_eclass_vtk_type[tree->eclass]);
    if (freturn <= 0) {
      return 0;
    }
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_level_kernel (t8_forest_t forest,
                                  t8_locidx_t ltree_id, t8_tree_t tree,
                                  t8_locidx_t element_index,
                                  FILE * vtufile, int *columns,
                                  void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    t8_eclass_scheme_t *ts;
    t8_element_t       *element;
    ts = forest->scheme->eclass_schemes[tree->eclass];
    element =
      (t8_element_t *) sc_array_index (&tree->elements, element_index);
    fprintf (vtufile, "%i ", t8_element_level (ts, element));
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_rank_kernel (t8_forest_t forest,
                                 t8_locidx_t ltree_id, t8_tree_t tree,
                                 t8_locidx_t element_index,
                                 FILE * vtufile, int *columns,
                                 void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    fprintf (vtufile, "%i ", forest->mpirank);
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_treeid_kernel (t8_forest_t forest,
                                   t8_locidx_t ltree_id, t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    fprintf (vtufile, "%lli ",
             (long long) ltree_id + forest->first_local_tree);
    *columns += 1;
  }
  return 1;
}

/* Iterate over all cells and write cell data to the file using
 * the cell_data_kernel as callback */
static int
t8_forest_vtk_write_cell_data (t8_forest_t forest, FILE * vtufile,
                               const char *dataname, const char *datatype,
                               const char *component_string,
                               int max_columns,
                               t8_forest_vtk_cell_data_kernel kernel)
{
  int                 freturn;
  int                 countcols;
  t8_tree_t           tree;
  t8_locidx_t         itree;
  t8_locidx_t         element_index;
  t8_locidx_t         num_local_trees;
  void               *data;

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = fprintf (vtufile, "        <DataArray type=\"%s\" "
                     "Name=\"%s\" %s format=\"ascii\">\n         ",
                     datatype, dataname, component_string);
  if (freturn <= 0) {
    return 0;
  }

  /* Call the kernel in initilization modus to possibly initialize the
   * data pointer */
  kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_INIT);
  /* We iterate over the trees and count each trees vertices,
   * we add this to the already counted vertices and write it to the file */
  /* TODO: replace with an element iterator */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, countcols = 0; itree < num_local_trees; itree++) {
    /* Get the tree that stores the elements */
    tree = t8_forest_get_tree (forest, itree);
    for (element_index = 0;
         element_index < (t8_locidx_t) tree->elements.elem_count;
         element_index++) {
      /* Execute the given callback on each element */
      if (!kernel (forest, itree, tree, element_index, vtufile, &countcols,
                   &data, T8_VTK_KERNEL_EXECUTE)) {
        /* call the kernel in clean-up modus */
        kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
        return 0;
      }
      /* After max_columns we break the line */
      if (!(countcols % max_columns)) {
        freturn = fprintf (vtufile, "\n         ");
        if (freturn <= 0) {
          /* call the kernel in clean-up modus */
          kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
          return 0;
        }
      }
    }                           /* element loop ends here */
    if (freturn <= 0) {
      /* call the kernel in clean-up modus */
      kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
      return 0;
    }
  }                             /* tree loop ends here */
  freturn = fprintf (vtufile, "\n        </DataArray>\n");
  if (freturn <= 0) {
    /* call the kernel in clean-up modus */
    kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
    return 0;
  }

  /* call the kernel in clean-up modus */
  kernel (NULL, 0, NULL, 0, NULL, NULL, &data, T8_VTK_KERNEL_CLEANUP);
  return 1;
}

/* Write the cell data to an open file stream.
 * Returns true on success and zero otherwise.
 * After completion the file will remain open, whether writing
 * cells was successful or not. */
static int
t8_forest_vtk_write_cells (t8_forest_t forest, FILE * vtufile, int write_id)
{
  int                 freturn;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (vtufile != NULL);

  /* Write the vertex coordinates */

  freturn = fprintf (vtufile, "      <Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "Position",
                                           T8_VTK_FLOAT_NAME,
                                           "NumberOfComponents=\"3\"",
                                           8,
                                           t8_forest_vtk_cells_vertices_kernel);
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  freturn = fprintf (vtufile, "      </Points>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done writing vertex coordinates */

  freturn = fprintf (vtufile, "      <Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "connectivity",
                                           T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_connectivity_kernel);
  if (freturn == 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the connectivity */

  /* Write the offsets, that is for each tree the index of the first entry
   * in the connectivity output that
   * does not refer to a vertex of the tree anymore.
   * For example if the trees are a square and a triangle, the offsets would
   * be 4 and 7, since indices 0,1,2,3 refer to the vertices of the square
   * and indices 4,5,6 to the indices of the triangle. */
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "offsets",
                                           T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_offset_kernel);
  if (freturn == 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the offsets */

  /* Write the element types. The type specifies the element class, thus
   * square/triangle/tet etc. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "types",
                                           "Int32", "", 8,
                                           t8_forest_vtk_cells_type_kernel);

  /* Done with writing the types */
  freturn = fprintf (vtufile, "      </Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  freturn = fprintf (vtufile, "      <CellData Scalars =\"%s%s\">\n",
                     "treeid,mpirank,level", (write_id ? "id" : ""));
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Write the tree ids. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "treeid",
                                           T8_VTK_GLOIDX, "", 8,
                                           t8_forest_vtk_cells_treeid_kernel);

  /* Done with writing the tree ids */

  /* Write the mpiranks. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "mpirank",
                                           "Int32", "", 8,
                                           t8_forest_vtk_cells_rank_kernel);

  /* Done with writing the mpiranks */

  /* Write the element refinement levels. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "level",
                                           "Int32", "", 8,
                                           t8_forest_vtk_cells_level_kernel);

  /* Done with writing the levels */
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

int
t8_forest_vtk_write_file (t8_forest_t forest, const char *fileprefix,
                          int write_id)
{
  FILE               *vtufile;
  t8_locidx_t         num_elements, num_points;
  char                vtufilename[BUFSIZ];
  int                 freturn;
  int                 write_treeid, write_rank, write_level;    /* TODO: these are set true currently
                                                                   make them optional */

  T8_ASSERT (forest != NULL);
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (fileprefix != NULL);
  T8_ASSERT (!write_id);        /* Write_id is not implemented */

  /* Currently we only support output in ascii format, not binary */
  T8_ASSERT (T8_VTK_ASCII == 1);

  /* process 0 creates the .pvtu file */
  if (forest->mpirank == 0) {
    write_treeid = write_rank = write_level = 1;
    if (t8_write_pvtu (fileprefix, forest->mpisize, write_treeid, write_rank,
                       write_level, write_id)) {
      t8_errorf ("Error when writing file %s.pvtu\n", fileprefix);
      goto t8_forest_vtk_failure;
    }
  }

  /* The local number of elements */
  num_elements = t8_forest_get_num_element (forest);
  /* The local number of points, counted with multiplicity */
  num_points = t8_forest_num_points (forest);

  /* The filename for this processes file */
  freturn =
    snprintf (vtufilename, BUFSIZ, "%s_%04d.vtu", fileprefix,
              forest->mpirank);
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
  freturn =
    fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
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
  freturn = fprintf (vtufile,
                     "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
                     (long long) num_points, (long long) num_elements);
  if (freturn <= 0) {
    goto t8_forest_vtk_failure;
  }

  /* write the cell data */
  if (!t8_forest_vtk_write_cells (forest, vtufile, write_id)) {
    /* Writing cells was not successful */
    goto t8_forest_vtk_failure;
  }

  freturn = fprintf (vtufile, "    </Piece>\n"
                     "  </UnstructuredGrid>\n" "</VTKFile>\n");
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

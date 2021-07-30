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
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_forest_types.h"
#include "t8_schemes/t8_default/t8_dpyramid_bits.h"
#if T8_WITH_VTK
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkPyramid.h>
#include <vtkWedge.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSmartPointer.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif
#include <t8.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_eclass.h>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

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
 *  CLEANUP - Called once after all elements. Used to cleanup any memory
 *            allocated during INIT.
 * \param [in] forest The forest.
 * \param [in] ltree_id   A local treeid.
 * \param [in] tree   The local tree of the forest with id \a ltree_id.
 * \param [in] element_index An index of an element inside \a tree.
 * \param [in] element  A pointer to the current element.
 * \param [in] ts       The eclass scheme of the current element.
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
typedef int         (*t8_forest_vtk_cell_data_kernel) (t8_forest_t forest,
                                                       t8_locidx_t ltree_id,
                                                       t8_tree_t tree,
                                                       t8_locidx_t
                                                       element_index,
                                                       t8_element_t * element,
                                                       t8_eclass_scheme_c *
                                                       ts, int is_ghost,
                                                       FILE * vtufile,
                                                       int *columns,
                                                       void **data,
                                                       T8_VTK_KERNEL_MODUS
                                                       modus);

void
t8_forest_write_vtk_via_API (t8_forest_t forest, const char *fileprefix)
{
#if T8_WITH_VTK
  /*Check assertions: forest and fileprefix are not NULL and forest is commited */
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);
  T8_ASSERT (fileprefix != NULL);

  long int            point_id = 0;     /* The id of the point in the points Object. */
  t8_locidx_t         ielement; /* The iterator over elements in a tree. */
  t8_locidx_t         itree, ivertex;
  double             *vertices;
  double              coordinates[3];
  int                 elem_id = 0;

/* Since we want to use different element types and a points Array and cellArray 
 * we have to declare these vtk objects. The cellArray stores the Elements.
 * The points and cellArray are needed to store the data we want to write in the Unstructured Grid. 
 */
  vtkNew < vtkPoints > points;
  vtkNew < vtkCellArray > cellArray;
  vtkNew < vtkHexahedron > hexa;
  vtkNew < vtkVertex > vertex;
  vtkNew < vtkLine > line;
  vtkNew < vtkQuad > quad;
  vtkNew < vtkTriangle > tri;
  vtkNew < vtkPyramid > pyramid;
  vtkNew < vtkWedge > prism;
  vtkNew < vtkTetra > tet;

  /* 
   * The cellTypes Array stores the element types as integers(see vtk doc).
   */
  int                *cellTypes =
    T8_ALLOC (int, t8_forest_get_local_num_elements (forest));

/* We iterate over all local trees*/
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
/* 
 * We get the vertices, the current tree, the scheme for this tree
 * and the number of elements in this tree. We need the vertices of
 * the tree to get the coordinates of the elements later. We need
 * the number of elements in this tree to iterate over all of them.
 */
    vertices = t8_forest_get_tree_vertices (forest, itree);
    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    t8_locidx_t         elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);

    /* We iterate over all elements in the tree */
    for (ielement = 0; ielement < elems_in_tree; ielement++) {

      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      vtkSmartPointer < vtkCell > pvtkCell = NULL;
      t8_element_shape_t  element_shape = scheme->t8_element_shape (element);
      int                 num_corners =
        scheme->t8_element_num_corners (element);

      /* depending on the element type we choose the correct vtk cell to insert points to */
      switch (element_shape) {
      case T8_ECLASS_VERTEX:
        pvtkCell = vertex;
        break;
      case T8_ECLASS_LINE:
        pvtkCell = line;
        break;
      case T8_ECLASS_QUAD:
        pvtkCell = quad;
        break;
      case T8_ECLASS_TRIANGLE:
        pvtkCell = tri;
        break;
      case T8_ECLASS_HEX:
        pvtkCell = hexa;
        break;
      case T8_ECLASS_TET:
        pvtkCell = tet;
        break;
      case T8_ECLASS_PRISM:
        pvtkCell = prism;
        break;
      case T8_ECLASS_PYRAMID:
        pvtkCell = pyramid;
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }

      /* For each element we iterate over all points */
      for (ivertex = 0; ivertex < num_corners; ivertex++, point_id++) {
        /* We take the element coordinates in vtk order */
        t8_forest_element_coordinate (forest, itree, element,
                                      vertices,
                                      t8_eclass_vtk_corner_number
                                      [element_shape]
                                      [ivertex], coordinates);

        /* Insert point in the points array */
        points->InsertNextPoint (coordinates[0], coordinates[1],
                                 coordinates[2]);
        /* Set the point ids to the vtk cell */
        pvtkCell->GetPointIds ()->SetId (ivertex, point_id);
      }
      /* We insert the next cell in the cell array */
      cellArray->InsertNextCell (pvtkCell);

      cellTypes[elem_id] = t8_eclass_vtk_type[element_shape];
      elem_id++;
    }
  }

  /* 
   * Write file: First we construct the unstructured Grid 
   * that will store the points and elements. It requires
   * information about the points(coordinates, stored in the points object)
   * and the cells(cellTypes and which points belong to this cell) 
   */

  vtkNew < vtkUnstructuredGrid > unstructuredGrid;
  unstructuredGrid->SetPoints (points);
  unstructuredGrid->SetCells (cellTypes, cellArray);
  /*
   * We define the filename used to write the pvtu and the vtu files.
   * The pwriterObj is of class XMLPUnstructuredGridWriter, the P in
   * XMLP is important: We want to write a vtu file for each process.
   * This class enables us to do exactly that. 
   */
  char                mpifilename[BUFSIZ];
  snprintf (mpifilename, BUFSIZ, "%s.pvtu", fileprefix);

  vtkSmartPointer < vtkXMLPUnstructuredGridWriter > pwriterObj =
    vtkSmartPointer < vtkXMLPUnstructuredGridWriter >::New ();
/*
 * Get/Set whether the appended data section is base64 encoded. 
 * If encoded, reading and writing will be slower, but the file 
 * will be fully valid XML and text-only. 
 * If not encoded, the XML specification will be violated, 
 * but reading and writing will be fast. The default is to do the encoding.
 * Documentation: https://vtk.org/doc/release/5.0/html/a02260.html#z3560_2
 * 
 * We set the filename of the pvtu file. The filenames of the vtu files
 * are given based on the name of the pvtu file and the process number.
 */
  pwriterObj->EncodeAppendedDataOff ();
  pwriterObj->SetFileName (mpifilename);

/*
 * Since we want to write multiple files, the processes 
 * have to communicate. Therefore, we define the communicator
 * vtk_comm and set it as the communicator. 
 * We have to set a controller for the pwriterObj, 
 * therefore we define the controller vtk_mpi_ctrl.
 */
#if T8_ENABLE_MPI
  vtkSmartPointer < vtkMPICommunicator > vtk_comm =
    vtkSmartPointer < vtkMPICommunicator >::New ();
  vtkMPICommunicatorOpaqueComm vtk_opaque_comm (&forest->mpicomm);
  vtk_comm->InitializeExternal (&vtk_opaque_comm);

  vtkSmartPointer < vtkMPIController > vtk_mpi_ctrl =
    vtkSmartPointer < vtkMPIController >::New ();
  vtk_mpi_ctrl->SetCommunicator (vtk_comm);

  pwriterObj->SetController (vtk_mpi_ctrl);
#endif
/*
 * We set the number of pieces as the number of mpi processes,
 * since we want to write a file for each process. We also
 * need to define a Start and EndPiece for the current
 * process. Then we can set the inputData for the writer:
 * We want to write the unstructured Grid, update the writer
 * and then write.
 * 
 * Note: We could write more than one file per process here, if desired.
 */
  pwriterObj->SetNumberOfPieces (forest->mpisize);
  pwriterObj->SetStartPiece (forest->mpirank);
  pwriterObj->SetEndPiece (forest->mpirank);
  pwriterObj->SetInputData (unstructuredGrid);
  pwriterObj->Update ();
  pwriterObj->Write ();
/* We have to free the allocated memory for the cellTypes Array. */
  T8_FREE (cellTypes);
#else
  t8_global_errorf
    ("Warning: t8code is not linked against vtk library. Vtk output will not be generated.\n");
  t8_global_productionf
    ("Consider calling 't8_forest_write_vtk' or 't8_forest_vtk_write_file' instead.\n");
#endif
}

static              t8_locidx_t
t8_forest_num_points (t8_forest_t forest, int count_ghosts)
{
  t8_locidx_t         itree, num_points, num_ghosts, ielem;
  t8_tree_t           tree;
  t8_eclass_t         ghost_class;
  size_t              num_elements;
  t8_dpyramid_t      *pyra;
  t8_element_array_t *ghost_elem;

  num_points = 0;
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* Get the tree that stores the elements */
    const t8_eclass_t   eclass = t8_forest_get_tree_class (forest, itree);
    t8_eclass_scheme_c *scheme = t8_forest_get_eclass_scheme (forest, eclass);
    num_elements = t8_forest_get_tree_num_elements (forest, itree);

    for (ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
      const t8_element_t *element =
        t8_forest_get_element_in_tree (forest, itree, ielem);
      const t8_element_shape_t shape = scheme->t8_element_shape (element);
      num_points += t8_eclass_num_vertices[shape];
    }
  }
  if (count_ghosts) {
    T8_ASSERT (forest->ghosts != NULL);
    /* We also count the points of the ghost cells */
    num_ghosts = t8_forest_ghost_num_trees (forest);
    for (itree = 0; itree < num_ghosts; itree++) {
      /* Get the element class of the ghost */
      ghost_class = t8_forest_ghost_get_tree_class (forest, itree);
      if (ghost_class == T8_ECLASS_PYRAMID) {
        ghost_elem = t8_forest_ghost_get_tree_elements (forest, itree);
        num_elements = t8_forest_ghost_tree_num_elements (forest, itree);
        for (ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
          pyra =
            (t8_dpyramid_t *) t8_element_array_index_locidx (ghost_elem,
                                                             ielem);
          num_points += t8_dpyramid_num_vertices (pyra);
        }
      }
      else {
        num_points += t8_eclass_num_vertices[ghost_class]
          * t8_forest_ghost_tree_num_elements (forest, itree);
      }
    }
  }
  return num_points;
}

static int
t8_forest_vtk_cells_vertices_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                     t8_tree_t tree,
                                     t8_locidx_t element_index,
                                     t8_element_t * element,
                                     t8_eclass_scheme_c * ts,
                                     int is_ghost,
                                     FILE * vtufile, int *columns,
                                     void **data, T8_VTK_KERNEL_MODUS modus)
{
  struct t8_forest_vtk_vertices_t
  {
    t8_locidx_t         ltreeid;        /* Store the last treeid with which the lernel was called.
                                           This is either a local tree id or a local ghost tree id */
    double              tree_vertices[T8_ECLASS_MAX_CORNERS * 3];       /* Stores the vertex coordinates of the tree */
  }                  *vertex_data;

#if 0
  /* if we eventually implement scaling the elements, activate this line */
  double              midpoint[3];
#endif
  double              element_coordinates[3];
  int                 num_tree_vertices, num_el_vertices, ivertex;
  int                 freturn;

  if (modus == T8_VTK_KERNEL_INIT) {
    /* We initialize the user data to store NULL as the current tree */
    *data = T8_ALLOC_ZERO (struct t8_forest_vtk_vertices_t, 1);
    vertex_data = (struct t8_forest_vtk_vertices_t *) *data;
    /* Set an invalid tree id as first tree id */
    vertex_data->ltreeid = -1;
    return 1;
  }
  else if (modus == T8_VTK_KERNEL_CLEANUP) {
    T8_FREE (*data);
    return 1;
  }
  num_tree_vertices = t8_eclass_num_vertices[ts->eclass];
  T8_ASSERT (modus == T8_VTK_KERNEL_EXECUTE);
  vertex_data = (struct t8_forest_vtk_vertices_t *) *data;
  if (ltree_id != vertex_data->ltreeid) {
    t8_cmesh_t          cmesh;
    t8_locidx_t         cmesh_local_id;
    double             *temp_vertices;

    /* The current tree is not the tree that we stored from
     * the last call to this function */
    vertex_data->ltreeid = ltree_id;
    /* get the coarse mesh tree */
    cmesh = t8_forest_get_cmesh (forest);
    /* Comput the cmesh local id of the tree */
    cmesh_local_id = t8_forest_ltreeid_to_cmesh_ltreeid (forest, ltree_id);
    /* Get the vertex coordinates of this tree */
    temp_vertices = ((double *) t8_cmesh_get_attribute (cmesh,
                                                        t8_get_package_id (),
                                                        0, cmesh_local_id));

    /* Copy the tree's vertex coordinates into the struct of the data pointer */
    memcpy (vertex_data->tree_vertices, temp_vertices, sizeof (*temp_vertices)
            * num_tree_vertices * 3);
  }

  /* TODO: be careful with pyramid class here.
   *       does this work too over tree->class or do we need something else?
   */

#if 0
  /* if we eventually implement scaling the elements, activate this line */
  t8_forest_element_centroid (forest, ltree_id, element,
                              vertex_data->tree_vertices, midpoint);
#endif
  num_el_vertices = t8_eclass_num_vertices[ts->t8_element_shape (element)];
  for (ivertex = 0; ivertex < num_el_vertices; ivertex++) {
    t8_forest_element_coordinate (forest, ltree_id, element,
                                  vertex_data->tree_vertices,
                                  t8_eclass_vtk_corner_number
                                  [ts->t8_element_shape (element)]
                                  [ivertex], element_coordinates);
#if 0
    /* if we eventually implement scaling the elements, activate this line */
    /* replace 0.9 with the scale factor
     * replace 0.1 with 1-scale_factor */
    t8_vec_ax (element_coordinates, 0.9);
    t8_vec_axpy (midpoint, element_coordinates, 0.1);
#endif
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
                                         t8_element_t * elements,
                                         t8_eclass_scheme_c * ts,
                                         int is_ghost,
                                         FILE * vtufile, int *columns,
                                         void **data,
                                         T8_VTK_KERNEL_MODUS modus)
{
  int                 ivertex, num_vertices;
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
  //SC_CHECK_ABORT (ts->eclass != T8_ECLASS_PYRAMID,
  //                "No vtk support for pyramids.");
  num_vertices = t8_eclass_num_vertices[ts->t8_element_shape (elements)];
  for (ivertex = 0; ivertex < num_vertices; ++ivertex, (*count_vertices)++) {
    freturn = fprintf (vtufile, " %ld", (long) *count_vertices);
    if (freturn <= 0) {
      return 0;
    }
  }
  *columns += t8_eclass_num_vertices[ts->t8_element_shape (elements)];
  return 1;
}

static int
t8_forest_vtk_cells_offset_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                   t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   t8_element_t * element,
                                   t8_eclass_scheme_c * ts,
                                   int is_ghost,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  long long          *offset;
  int                 freturn;
  int                 num_vertices;

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
  //SC_CHECK_ABORT (ts->eclass != T8_ECLASS_PYRAMID,
  //                "Pyramids not supported in vtk");
  num_vertices = t8_eclass_num_vertices[ts->t8_element_shape (element)];
  *offset += num_vertices;
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
                                 t8_element_t * element,
                                 t8_eclass_scheme_c * ts,
                                 int is_ghost,
                                 FILE * vtufile, int *columns,
                                 void **data, T8_VTK_KERNEL_MODUS modus)
{
  int                 freturn;
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    /* print the vtk type of the element */
    freturn =
      fprintf (vtufile, " %d",
               t8_eclass_vtk_type[ts->t8_element_shape (element)]);
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
                                  t8_element_t * element,
                                  t8_eclass_scheme_c * ts,
                                  int is_ghost,
                                  FILE * vtufile, int *columns,
                                  void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    fprintf (vtufile, "%i ", ts->t8_element_level (element));
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_rank_kernel (t8_forest_t forest,
                                 t8_locidx_t ltree_id, t8_tree_t tree,
                                 t8_locidx_t element_index,
                                 t8_element_t * element,
                                 t8_eclass_scheme_c * ts,
                                 int is_ghost,
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
                                   t8_element_t * element,
                                   t8_eclass_scheme_c * ts,
                                   int is_ghost,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    long long           tree_id;
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
t8_forest_vtk_cells_elementid_kernel (t8_forest_t forest,
                                      t8_locidx_t ltree_id, t8_tree_t tree,
                                      t8_locidx_t element_index,
                                      t8_element_t * element,
                                      t8_eclass_scheme_c * ts,
                                      int is_ghost,
                                      FILE * vtufile, int *columns,
                                      void **data, T8_VTK_KERNEL_MODUS modus)
{
  if (modus == T8_VTK_KERNEL_EXECUTE) {
    if (!is_ghost) {
      fprintf (vtufile, "%lli ", element_index + tree->elements_offset +
               (long long) t8_forest_get_first_local_element_id (forest));
    }
    else {
      fprintf (vtufile, "%lli ", (long long) -1);
    }
    *columns += 1;
  }
  return 1;
}

static int
t8_forest_vtk_cells_scalar_kernel (t8_forest_t forest,
                                   t8_locidx_t ltree_id, t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   t8_element_t * element,
                                   t8_eclass_scheme_c * ts,
                                   int is_ghost,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  double              element_value = 0;
  t8_locidx_t         scalar_index;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    /* For local elements access the data array, for ghosts, write 0 */
    if (!is_ghost) {
      scalar_index = t8_forest_get_tree_element_offset (forest, ltree_id)
        + element_index;
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
t8_forest_vtk_cells_vector_kernel (t8_forest_t forest,
                                   t8_locidx_t ltree_id, t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   t8_element_t * element,
                                   t8_eclass_scheme_c * ts,
                                   int is_ghost,
                                   FILE * vtufile, int *columns,
                                   void **data, T8_VTK_KERNEL_MODUS modus)
{
  double             *element_values, null_vec[3] = { 0, 0, 0 };
  int                 dim, idim;
  t8_locidx_t         tree_offset;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    dim = 3;
    T8_ASSERT (forest->dimension <= 3);
    /* For local elements access the data array, for ghosts, write 0 */
    if (!is_ghost) {
      tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);
      /* Get a pointer to the start of the element's vector data */
      element_values =
        ((double *) *data) + (tree_offset + element_index) * dim;
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
t8_forest_vtk_vertices_scalar_kernel (t8_forest_t forest,
                                      t8_locidx_t ltree_id, t8_tree_t tree,
                                      t8_locidx_t element_index,
                                      t8_element_t * element,
                                      t8_eclass_scheme_c * ts,
                                      int is_ghost,
                                      FILE * vtufile, int *columns,
                                      void **data, T8_VTK_KERNEL_MODUS modus)
{
  double              element_value = 0;
  int                 num_vertex, ivertex;
  t8_locidx_t         scalar_index;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    num_vertex = ts->t8_element_num_corners (element);

    for (ivertex = 0; ivertex < num_vertex; ivertex++) {
      /* For local elements access the data array, for ghosts, write 0 */
      if (!is_ghost) {
        scalar_index = t8_forest_get_tree_element_offset (forest, ltree_id)
          + element_index;
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
t8_forest_vtk_vertices_vector_kernel (t8_forest_t forest,
                                      t8_locidx_t ltree_id, t8_tree_t tree,
                                      t8_locidx_t element_index,
                                      t8_element_t * element,
                                      t8_eclass_scheme_c * ts,
                                      int is_ghost,
                                      FILE * vtufile, int *columns,
                                      void **data, T8_VTK_KERNEL_MODUS modus)
{
  double             *element_values, null_vec[3] = { 0, 0, 0 };
  int                 dim, idim;
  int                 num_vertex, ivertex;
  t8_locidx_t         tree_offset;

  if (modus == T8_VTK_KERNEL_EXECUTE) {
    num_vertex = ts->t8_element_num_corners (element);
    for (ivertex = 0; ivertex < num_vertex; ivertex++) {
      dim = 3;
      T8_ASSERT (forest->dimension <= 3);
      /* For local elements access the data array, for ghosts, write 0 */
      if (!is_ghost) {
        tree_offset = t8_forest_get_tree_element_offset (forest, ltree_id);
        /* Get a pointer to the start of the element's vector data */
        element_values =
          ((double *) *data) + (tree_offset + element_index) * dim;
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
t8_forest_vtk_write_cell_data (t8_forest_t forest, FILE * vtufile,
                               const char *dataname,
                               const char *datatype,
                               const char *component_string,
                               int max_columns,
                               t8_forest_vtk_cell_data_kernel kernel,
                               int write_ghosts, void *udata)
{
  int                 freturn;
  int                 countcols;
  t8_tree_t           tree;
  t8_locidx_t         itree, ighost;
  t8_locidx_t         element_index, elems_in_tree;
  t8_locidx_t         num_local_trees, num_ghost_trees;
  t8_element_t       *element;
  t8_eclass_scheme_c *ts;
  void               *data = NULL;

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = fprintf (vtufile, "        <DataArray type=\"%s\" "
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

  /* Call the kernel in initilization modus to possibly initialize the
   * data pointer */
  kernel (NULL, 0, NULL, 0, NULL, NULL, 0, NULL, NULL, &data,
          T8_VTK_KERNEL_INIT);
  /* We iterate over the trees and count each trees vertices,
   * we add this to the already counted vertices and write it to the file */
  /* TODO: replace with an element iterator */
  num_local_trees = t8_forest_get_num_local_trees (forest);
  for (itree = 0, countcols = 0; itree < num_local_trees; itree++) {
    /* Get the tree that stores the elements */
    tree = t8_forest_get_tree (forest, itree);
    /* Get the eclass scheme of the tree */
    ts = t8_forest_get_eclass_scheme (forest,
                                      t8_forest_get_tree_class (forest,
                                                                itree));
    elems_in_tree =
      (t8_locidx_t) t8_element_array_get_count (&tree->elements);
    for (element_index = 0; element_index < elems_in_tree; element_index++) {
      /* Get a pointer to the element */
      element =
        t8_forest_get_element (forest, tree->elements_offset + element_index,
                               NULL);
      T8_ASSERT (element != NULL);
      /* Execute the given callback on each element */
      if (!kernel
          (forest, itree, tree, element_index, element, ts, 0, vtufile,
           &countcols, &data, T8_VTK_KERNEL_EXECUTE)) {
        /* call the kernel in clean-up modus */
        kernel (NULL, 0, NULL, 0, NULL, NULL, 0, NULL, NULL, &data,
                T8_VTK_KERNEL_CLEANUP);
        return 0;
      }
      /* After max_columns we break the line */
      if (!(countcols % max_columns)) {
        freturn = fprintf (vtufile, "\n         ");
        if (freturn <= 0) {
          /* call the kernel in clean-up modus */
          kernel (NULL, 0, NULL, 0, NULL, NULL, 0, NULL, NULL, &data,
                  T8_VTK_KERNEL_CLEANUP);
          return 0;
        }
      }
    }                           /* element loop ends here */
    if (freturn <= 0) {
      /* call the kernel in clean-up modus */
      kernel (NULL, 0, NULL, 0, NULL, NULL, 0, NULL, NULL, &data,
              T8_VTK_KERNEL_CLEANUP);
      return 0;
    }
  }                             /* tree loop ends here */

  if (write_ghosts) {
    t8_locidx_t         num_ghosts_in_tree;
    /* Iterate over the ghost elements */
    /* TODO: replace with an element iterator */
    num_ghost_trees = t8_forest_ghost_num_trees (forest);
    for (ighost = 0; ighost < num_ghost_trees; ighost++) {
      /* Get the eclass scheme of the ghost tree */
      ts =
        t8_forest_get_eclass_scheme (forest,
                                     t8_forest_ghost_get_tree_class (forest,
                                                                     ighost));
      /* The number of ghosts in this tree */
      num_ghosts_in_tree = t8_forest_ghost_tree_num_elements (forest, ighost);
      for (element_index = 0;
           element_index < num_ghosts_in_tree; element_index++) {
        /* Get a pointer to the element */
        element = t8_forest_ghost_get_element (forest, ighost, element_index);
        /* Execute the given callback on each element */
        if (!kernel
            (forest, ighost + num_local_trees, NULL, element_index, element,
             ts, 1, vtufile, &countcols, &data, T8_VTK_KERNEL_EXECUTE)) {
          /* call the kernel in clean-up modus */
          kernel (NULL, 0, NULL, 0, NULL, NULL, 1, NULL, NULL, &data,
                  T8_VTK_KERNEL_CLEANUP);
          return 0;
        }
        /* After max_columns we break the line */
        if (!(countcols % max_columns)) {
          freturn = fprintf (vtufile, "\n         ");
          if (freturn <= 0) {
            /* call the kernel in clean-up modus */
            kernel (NULL, 0, NULL, 0, NULL, NULL, 1, NULL, NULL, &data,
                    T8_VTK_KERNEL_CLEANUP);
            return 0;
          }
        }
      }                         /* element loop ends here */
      if (freturn <= 0) {
        /* call the kernel in clean-up modus */
        kernel (NULL, 0, NULL, 0, NULL, NULL, 1, NULL, NULL, &data,
                T8_VTK_KERNEL_CLEANUP);
        return 0;
      }
    }                           /* ghost loop ends here */
  }                             /* write_ghosts ends here */
  /* call the kernel in clean-up modus */
  kernel (NULL, 0, NULL, 0, NULL, NULL, 0, NULL, NULL, &data,
          T8_VTK_KERNEL_CLEANUP);
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
t8_forest_vtk_write_cells (t8_forest_t forest, FILE * vtufile,
                           int write_treeid,
                           int write_mpirank,
                           int write_level, int write_element_id,
                           int write_ghosts, int num_data,
                           t8_vtk_data_field_t * data)
{
  int                 freturn;
  int                 idata;

  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (vtufile != NULL);

  freturn = fprintf (vtufile, "      <Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  /* Write the connectivity information.
   * Thus for each tree we write the indices of its corner vertices. */
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "connectivity",
                                           T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_connectivity_kernel,
                                           write_ghosts, NULL);
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
  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "offsets",
                                           T8_VTK_LOCIDX, "", 8,
                                           t8_forest_vtk_cells_offset_kernel,
                                           write_ghosts, NULL);
  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the offsets */

  /* Write the element types. The type specifies the element class, thus
   * square/triangle/tet etc. */

  freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "types",
                                           "Int32", "", 8,
                                           t8_forest_vtk_cells_type_kernel,
                                           write_ghosts, NULL);

  if (!freturn) {
    goto t8_forest_vtk_cell_failure;
  }
  /* Done with writing the types */
  freturn = fprintf (vtufile, "      </Cells>\n");
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  freturn = fprintf (vtufile, "      <CellData Scalars =\"%s%s\">\n",
                     "treeid,mpirank,level", (write_element_id ? "id" : ""));
  if (freturn <= 0) {
    goto t8_forest_vtk_cell_failure;
  }

  if (write_treeid) {
    /* Write the tree ids. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "treeid",
                                             T8_VTK_GLOIDX, "", 8,
                                             t8_forest_vtk_cells_treeid_kernel,
                                             write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the tree ids */
  }
  if (write_mpirank) {
    /* Write the mpiranks. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "mpirank",
                                             "Int32", "", 8,
                                             t8_forest_vtk_cells_rank_kernel,
                                             write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the mpiranks */
  }
  if (write_level) {
    /* Write the element refinement levels. */

    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "level",
                                             "Int32", "", 8,
                                             t8_forest_vtk_cells_level_kernel,
                                             write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }
    /* Done with writing the levels */
  }

  if (write_element_id) {
    /* Write the element ids. */
    const char         *datatype;

    /* Use 32 bit ints if the global element count fits, 64 bit otherwise. */
    datatype = forest->global_num_elements > T8_LOCIDX_MAX ? T8_VTK_GLOIDX :
      T8_VTK_LOCIDX;
    freturn = t8_forest_vtk_write_cell_data (forest, vtufile, "element_id",
                                             datatype, "", 8,
                                             t8_forest_vtk_cells_elementid_kernel,
                                             write_ghosts, NULL);
    if (!freturn) {
      goto t8_forest_vtk_cell_failure;
    }

    /* Done with writing the element ids */
  }
  /* Write the user defined data fields per element */
  for (idata = 0; idata < num_data; idata++) {
    if (data[idata].type == T8_VTK_SCALAR) {
      freturn =
        t8_forest_vtk_write_cell_data (forest, vtufile,
                                       data[idata].description,
                                       T8_VTK_FLOAT_NAME, "", 8,
                                       t8_forest_vtk_cells_scalar_kernel,
                                       write_ghosts, data[idata].data);
    }
    else {
      char                component_string[BUFSIZ];
      T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
      snprintf (component_string, BUFSIZ, "NumberOfComponents=\"3\"");
      freturn =
        t8_forest_vtk_write_cell_data (forest, vtufile,
                                       data[idata].description,
                                       T8_VTK_FLOAT_NAME,
                                       component_string,
                                       8 * forest->dimension,
                                       t8_forest_vtk_cells_vector_kernel,
                                       write_ghosts, data[idata].data);
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
t8_forest_vtk_write_points (t8_forest_t forest, FILE * vtufile,
                            int write_ghosts,
                            int num_data, t8_vtk_data_field_t * data)
{
  int                 freturn;
  int                 sreturn;
  int                 idata;
  char                description[BUFSIZ];

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
                                           t8_forest_vtk_cells_vertices_kernel,
                                           write_ghosts, NULL);
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
        sreturn =
          snprintf (description, BUFSIZ, "%s_%s", data[idata].description,
                    "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          t8_debugf
            ("Warning: Truncated vtk point data description to '%s'\n",
             description);
        }
        freturn =
          t8_forest_vtk_write_cell_data (forest, vtufile, description,
                                         T8_VTK_FLOAT_NAME, "", 8,
                                         t8_forest_vtk_vertices_scalar_kernel,
                                         write_ghosts, data[idata].data);
      }
      else {
        char                component_string[BUFSIZ];
        T8_ASSERT (data[idata].type == T8_VTK_VECTOR);
        snprintf (component_string, BUFSIZ, "NumberOfComponents=\"3\"");
        /* Write the description string. */
        sreturn =
          snprintf (description, BUFSIZ, "%s_%s", data[idata].description,
                    "points");

        if (sreturn >= BUFSIZ) {
          /* The output was truncated */
          /* Note: gcc >= 7.1 prints a warning if we 
           * do not check the return value of snprintf. */
          t8_debugf
            ("Warning: Truncated vtk point data description to '%s'\n",
             description);
        }

        freturn =
          t8_forest_vtk_write_cell_data (forest, vtufile, description,
                                         T8_VTK_FLOAT_NAME, component_string,
                                         8 * forest->dimension,
                                         t8_forest_vtk_vertices_vector_kernel,
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
t8_forest_vtk_write_file (t8_forest_t forest, const char *fileprefix,
                          int write_treeid,
                          int write_mpirank,
                          int write_level, int write_element_id,
                          int write_ghosts,
                          int num_data, t8_vtk_data_field_t * data)
{
  FILE               *vtufile = NULL;
  t8_locidx_t         num_elements, num_points;
  char                vtufilename[BUFSIZ];
  int                 freturn;

  T8_ASSERT (forest != NULL);
  T8_ASSERT (t8_forest_is_committed (forest));
  T8_ASSERT (fileprefix != NULL);
  if (forest->ghosts == NULL || forest->ghosts->num_ghosts_elements == 0) {
    /* Never write ghost elements if there aren't any */
    write_ghosts = 0;
  }
  T8_ASSERT (forest->ghosts != NULL || !write_ghosts);

  /* Currently we only support output in ascii format, not binary */
  T8_ASSERT (T8_VTK_ASCII == 1);

  /* process 0 creates the .pvtu file */
  if (forest->mpirank == 0) {
    if (t8_write_pvtu
        (fileprefix, forest->mpisize, write_treeid, write_mpirank,
         write_level, write_element_id, num_data, data)) {
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
  /* write the point data */
  if (!t8_forest_vtk_write_points
      (forest, vtufile, write_ghosts, num_data, data)) {
    /* writings points was not succesful */
    goto t8_forest_vtk_failure;
  }
  /* write the cell data */
  if (!t8_forest_vtk_write_cells
      (forest, vtufile, write_treeid, write_mpirank, write_level,
       write_element_id, write_ghosts, num_data, data)) {
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

T8_EXTERN_C_END ();

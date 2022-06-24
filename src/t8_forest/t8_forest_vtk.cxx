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

#include <t8_forest/t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"
#include "t8_forest_types.h"
#if T8_WITH_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkVertex.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkPyramid.h>
#include <vtkWedge.h>
#include <vtkQuadraticEdge.h>
#include <vtkQuadraticTriangle.h>
#include <vtkQuadraticQuad.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkQuadraticWedge.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkTypeInt64Array.h>
#if T8_ENABLE_MPI
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#endif
#endif
#include <t8.h>
#include <t8_forest.h>

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
                                                       t8_element_t *element,
                                                       t8_eclass_scheme_c *ts,
                                                       int is_ghost,
                                                       FILE *vtufile,
                                                       int *columns,
                                                       void **data,
                                                       T8_VTK_KERNEL_MODUS
                                                       modus);

#if T8_WITH_VTK
/* lookup table for number of nodes for curved eclasses. */
const int           t8_curved_eclass_num_nodes[T8_ECLASS_COUNT] =
  { 1, 3, 8, 6, 20, 10, 15, 13 };

/* lookup table for vtk types of curved elements */
const int           t8_curved_eclass_vtk_type[T8_ECLASS_COUNT] =
  { 1, 21, 23, 22, 25, 24, 26, 27 };
#endif

/* 
 * depending on whether we want to write curved or non-curved elements
 * we need the right number of points, so we choose the right lookup table
 */
#if T8_WITH_VTK
static int
t8_get_number_of_vtk_nodes (t8_element_shape_t eclass, int curved_flag)
{
  /* use the lookup table of the eclasses. */
  if (curved_flag) {
    return t8_curved_eclass_num_nodes[eclass];
  }
  return t8_eclass_num_vertices[eclass];
}
#endif

/* If we want to write curved elements, we need to calculate 
 * the reference coordinates. For the vertices(end points)
 * of the elements, we can use t8_element_vertex_reference_coords 
 * to get them. But for curved elements, we also need nodes at the 
 * middle points of lines of elements. We get those coordinates by 
 * adding the vertices and multiplying by 0.5. To get the 
 * correct node, we use e.g. (vertex - 3) % 4, for each 
 * element there is a correct order, therefore we have those
 * formulas. For more information look into the vtk documentation.
 * TODO: Add Pyramids when they are merged into the dev branch.
 * */
#if T8_WITH_VTK
static void
t8_curved_element_get_reference_node_coords (const t8_element_t *elem,
                                             t8_element_shape_t eclass,
                                             t8_eclass_scheme_c *scheme,
                                             int vertex, double *coords)
{
  double              vertex_coords[3] = { 0, 0, 0 };
  int                 i;
  int                 j;

  switch (eclass) {
  case T8_ECLASS_VERTEX:
    scheme->t8_element_vertex_reference_coords (elem,
                                                t8_eclass_vtk_corner_number
                                                [eclass][vertex], coords);
    break;
  case T8_ECLASS_LINE:
    if (vertex < 2) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex - 1],
                                                  vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex - 2],
                                                  coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_QUAD:
    if (vertex < 4) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 4) % 4];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }

    break;
  case T8_ECLASS_TRIANGLE:
    if (0 <= vertex && vertex <= 2) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else {
      i = (vertex - 3) % 3;
      j = (vertex - 2) % 3;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_HEX:
    if (vertex < 8) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (7 < vertex && vertex < 12) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 8) % 4];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 7) % 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else if (11 < vertex && vertex < 16) {
      i = t8_eclass_vtk_corner_number[eclass][((vertex - 8) % 4) + 4];
      j = t8_eclass_vtk_corner_number[eclass][((vertex - 7) % 4) + 4];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 16];
      j = i + 4;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }

    break;
  case T8_ECLASS_TET:
    if (vertex < 4) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (3 < vertex && vertex < 7) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 4) % 3];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 7];
      j = 3;
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  case T8_ECLASS_PRISM:
    if (vertex < 6) {
      scheme->t8_element_vertex_reference_coords (elem,
                                                  t8_eclass_vtk_corner_number
                                                  [eclass][vertex], coords);
    }
    else if (5 < vertex && vertex < 9) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex - 3) % 3];
      j = t8_eclass_vtk_corner_number[eclass][(vertex - 2) % 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else if (8 < vertex && vertex < 12) {
      i = t8_eclass_vtk_corner_number[eclass][(vertex % 3) + 3];
      j = t8_eclass_vtk_corner_number[eclass][((vertex + 1) % 3) + 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    else {
      i = t8_eclass_vtk_corner_number[eclass][vertex % 12];
      j = t8_eclass_vtk_corner_number[eclass][(vertex % 12) + 3];
      scheme->t8_element_vertex_reference_coords (elem, i, vertex_coords);
      scheme->t8_element_vertex_reference_coords (elem, j, coords);
      /* Compute the average of those coordinates */
      t8_vec_axpy (vertex_coords, coords, 1);
      t8_vec_ax (coords, 0.5);
    }
    break;
  default:
    scheme->t8_element_vertex_reference_coords (elem,
                                                t8_eclass_vtk_corner_number
                                                [eclass][vertex], coords);
    break;
  }
}
#endif

int
t8_forest_vtk_write_file_via_API (t8_forest_t forest, const char *fileprefix,
                                  int write_treeid,
                                  int write_mpirank,
                                  int write_level,
                                  int write_element_id,
                                  int curved_flag,
                                  int num_data, t8_vtk_data_field_t *data)
{
#if T8_WITH_VTK
  /*Check assertions: forest and fileprefix are not NULL and forest is commited */
  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->rc.refcount > 0);
  T8_ASSERT (forest->committed);
  T8_ASSERT (fileprefix != NULL);

  long int            point_id = 0;     /* The id of the point in the points Object. */
  t8_locidx_t         ielement; /* The iterator over elements in a tree. */
  t8_gloidx_t         itree, ivertex;
  double              coordinates[3];
  double              vertex_coords[3] = { 0, 0, 0 };
  int                 elem_id = 0;
  t8_locidx_t         num_elements;
  int                 freturn = 0;
  t8_gloidx_t         gtreeid;
  t8_cmesh_t          cmesh;
  int                 num_node;

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
  vtkNew < vtkWedge > prism;
  vtkNew < vtkTetra > tet;
  vtkNew < vtkPyramid > pyra;

  vtkNew < vtkQuadraticEdge > quadraticedge;
  vtkNew < vtkQuadraticTriangle > quadratictri;
  vtkNew < vtkQuadraticQuad > quadraticquad;
  vtkNew < vtkQuadraticTetra > quadratictet;
  vtkNew < vtkQuadraticHexahedron > quadratichexa;
  vtkNew < vtkQuadraticWedge > quadraticprism;

  /* 
   * The cellTypes Array stores the element types as integers(see vtk doc).
   */
  num_elements = t8_forest_get_local_num_elements (forest);
  int                *cellTypes = T8_ALLOC (int, num_elements);

  /*
   * We need the vertex coords array to be of the 
   * correct dim. Since it is always the same
   * in one mesh, we take the dim of one element.
   * We add 1 if we look at a vertext(dim=0) because 
   * an array of size 0 is not allowed. 
   * Then we allocate memory, because we do not know
   * beforehand how many entries the array needs.
   */

  /*
   * We have to define the vtkTypeInt64Array that hold 
   * metadata if wanted. 
   */

  t8_vtk_gloidx_array_type_t *vtk_treeid = t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_mpirank =
    t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_level = t8_vtk_gloidx_array_type_t::New ();
  t8_vtk_gloidx_array_type_t *vtk_element_id =
    t8_vtk_gloidx_array_type_t::New ();

/*
 * We need the dataArray for writing double valued user defined data in the vtu files.
 * We want to write num_data many timesteps/arrays.
 * We need num_data many vtkDoubleArrays, so we need to allocate storage.
 * Later we call the constructor with: dataArrays[idata]=vtkDoubleArray::New()
 */
  vtkDoubleArray    **dataArrays;
  dataArrays = T8_ALLOC (vtkDoubleArray *, num_data);

  cmesh = t8_forest_get_cmesh (forest);
/* We iterate over all local trees*/
  for (itree = 0; itree < t8_forest_get_num_local_trees (forest); itree++) {
/* 
 * We get the current tree, the scheme for this tree
 * and the number of elements in this tree. We need the vertices of
 * the tree to get the coordinates of the elements later. We need
 * the number of elements in this tree to iterate over all of them.
 */

    t8_eclass_scheme_c *scheme =
      t8_forest_get_eclass_scheme (forest, t8_forest_get_tree_class (forest,
                                                                     itree));
    t8_locidx_t         elems_in_tree =
      t8_forest_get_tree_num_elements (forest, itree);
    t8_locidx_t         offset =
      t8_forest_get_tree_element_offset (forest, itree);
    /* We iterate over all elements in the tree */
    /* Compute the global tree id */
    gtreeid = t8_forest_global_tree_id (forest, itree);
    for (ielement = 0; ielement < elems_in_tree; ielement++) {
      t8_element_t       *element =
        t8_forest_get_element_in_tree (forest, itree, ielement);
      T8_ASSERT (element != NULL);
      vtkSmartPointer < vtkCell > pvtkCell = NULL;
      t8_element_shape_t  element_shape = scheme->t8_element_shape (element);
      num_node = t8_get_number_of_vtk_nodes (element_shape, curved_flag);
      /* depending on the element type we choose the correct vtk cell to insert points to */
      if (curved_flag == 0) {
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
          pvtkCell = pyra;
          break;
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }
      else {                    /* curved_flag != 0 */
        switch (element_shape) {
        case T8_ECLASS_VERTEX:
          pvtkCell = vertex;
          break;
        case T8_ECLASS_LINE:
          pvtkCell = quadraticedge;
          break;
        case T8_ECLASS_QUAD:
          pvtkCell = quadraticquad;
          break;
        case T8_ECLASS_TRIANGLE:
          pvtkCell = quadratictri;
          break;
        case T8_ECLASS_HEX:
          pvtkCell = quadratichexa;
          break;
        case T8_ECLASS_TET:
          pvtkCell = quadratictet;
          break;
        case T8_ECLASS_PRISM:
          pvtkCell = quadraticprism;
          break;
        case T8_ECLASS_PYRAMID:
          SC_CHECK_ABORT (element_shape != T8_ECLASS_PYRAMID,
                          "Quadratic Pyramids are not supported in vtk output");
        default:
          SC_ABORT_NOT_REACHED ();
        }
      }

      /* For each element we iterate over all points */
      for (ivertex = 0; ivertex < num_node; ivertex++, point_id++) {
        /* Compute the vertex coordinates inside [0,1]^dim reference cube. */
        if (curved_flag) {
          t8_curved_element_get_reference_node_coords (element, element_shape,
                                                       scheme, ivertex,
                                                       vertex_coords);
        }
        else {
          scheme->t8_element_vertex_reference_coords (element,
                                                      t8_eclass_vtk_corner_number
                                                      [element_shape]
                                                      [ivertex],
                                                      vertex_coords);
        }

        /* Evaluate the geometry */
        t8_geometry_evaluate (cmesh, gtreeid, vertex_coords, coordinates);

        /* Insert point in the points array */
        points->InsertNextPoint (coordinates[0], coordinates[1],
                                 coordinates[2]);
        /* Set the point ids to the vtk cell */
        pvtkCell->GetPointIds ()->SetId (ivertex, point_id);
      }                         /* end loop over all vertices of the element */

      /* We insert the next cell in the cell array */
      cellArray->InsertNextCell (pvtkCell);
      /*
       * Write current cell Type in the cell Types array at the elem_id index.
       * Depending on the values of the binary inputs write_treeid, 
       * write_mpirank and write_element_id we also fill the corresponding
       * arrays with the data we want(treeid,mpirank,element_id).
       * To get the element id, we have to add the local id in the tree 
       * plus theo
       */

      /* *INDENT-OFF* */
      if(curved_flag==0){
        cellTypes[elem_id] = t8_eclass_vtk_type[element_shape];
      }
      else{
        cellTypes[elem_id] = t8_curved_eclass_vtk_type[element_shape];
      }
      if (write_treeid == 1) {
        vtk_treeid->InsertNextValue (itree);
      }
      if (write_mpirank == 1) {
        vtk_mpirank->InsertNextValue (forest->mpirank);
      }
      if (write_level == 1) {
        vtk_level->InsertNextValue (scheme->t8_element_level (element));
      }
      if (write_element_id == 1) {
        vtk_element_id->InsertNextValue (elem_id + offset +
                                         t8_forest_get_first_local_element_id
                                         (forest));
      }
      /* *INDENT-ON* */
      elem_id++;
    }                           /* end of loop over elements */
  }                             /* end of loop over local trees */

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
  /* *INDENT-OFF* */
  if (write_treeid) {
    vtk_treeid->SetName ("treeid");
    unstructuredGrid->GetCellData()->AddArray(vtk_treeid);
  }
  if (write_mpirank) {
    vtk_mpirank->SetName ("mpirank");
    unstructuredGrid->GetCellData()->AddArray(vtk_mpirank);
  }
  if (write_level) {
    vtk_level->SetName ("level");
    unstructuredGrid->GetCellData()->AddArray(vtk_level);
  }
  if (write_element_id) {
    vtk_element_id->SetName ("element_id");
    unstructuredGrid->GetCellData()->AddArray(vtk_element_id);
  }
  /* *INDENT-ON* */

/* Write the user defined data fields. 
 * For that we iterate over the idata, set the name, the array
 * and then give this data to the unstructured Grid Object.
 * We differentiate between scalar and vector data.
 */
  for (int idata = 0; idata < num_data; idata++) {
    dataArrays[idata] = vtkDoubleArray::New ();
    if (data[idata].type == T8_VTK_SCALAR) {
      dataArrays[idata]->SetName (data[idata].description);     /* Set the name of the array */
      dataArrays[idata]->SetVoidArray (data[idata].data, num_elements, 1);      /* We write the data in the array from the input array */
      unstructuredGrid->GetCellData ()->AddArray (dataArrays[idata]);   /* We add the array to the cell data object */
    }
    else {
      dataArrays[idata]->SetName (data[idata].description);     /* Set the name of the array */
      dataArrays[idata]->SetNumberOfTuples (num_elements);      /* We want number of tuples=number of elements */
      dataArrays[idata]->SetNumberOfComponents (3);     /* Each tuples has 3 values */
      dataArrays[idata]->SetVoidArray (data[idata].data, num_elements * 3, 1);  /*  */
      unstructuredGrid->GetCellData ()->SetVectors (dataArrays[idata]); /*  */
    }
  }

/* We set the input data and write the vtu files. */
  pwriterObj->SetInputData (unstructuredGrid);
  pwriterObj->Update ();
  if (pwriterObj->Write ()) {
    /* Writing was successful */
    freturn = 1;
  }
  else {
    t8_errorf ("Error when writing vtk file.\n");
  }

/* We have to free the allocated memory for the cellTypes Array and the other arrays we allocated memory for. */

  vtk_treeid->Delete ();
  vtk_mpirank->Delete ();
  vtk_level->Delete ();
  vtk_element_id->Delete ();
  for (int idata = 0; idata < num_data; idata++) {
    dataArrays[idata]->Delete ();
  }

  T8_FREE (cellTypes);
  T8_FREE (dataArrays);
/* Return whether writing was successful */
  return freturn;

#else
  t8_global_errorf
    ("Warning: t8code is not linked against vtk library. Vtk output will not be generated.\n");
  t8_global_productionf
    ("Consider calling 't8_forest_write_vtk' or 't8_forest_vtk_write_file' instead.\n");
  return 0;
#endif
}

static t8_locidx_t
t8_forest_num_points (t8_forest_t forest, int count_ghosts)
{
  t8_locidx_t         itree, num_points, num_ghosts, ielem;
  t8_tree_t           tree;
  t8_eclass_t         ghost_class;
  size_t              num_elements;
  t8_element_t       *elem;
  t8_element_array_t *ghost_elem;
  t8_eclass_scheme   *tscheme;

  num_points = 0;
  for (itree = 0; itree < (t8_locidx_t) forest->trees->elem_count; itree++) {
    /* Get the tree that stores the elements */
    tree = (t8_tree_t) t8_sc_array_index_topidx (forest->trees, itree);
    /* Get the scheme of the current tree */
    tscheme = t8_forest_get_eclass_scheme (forest, tree->eclass);
    num_elements = t8_element_array_get_count (&tree->elements);
    for (ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
      elem = t8_element_array_index_locidx (&tree->elements, ielem);
      num_points += tscheme->t8_element_num_corners (elem);
    }

  }
  if (count_ghosts) {
    T8_ASSERT (forest->ghosts != NULL);
    /* We also count the points of the ghost cells */
    num_ghosts = t8_forest_ghost_num_trees (forest);
    for (itree = 0; itree < num_ghosts; itree++) {
      /* Get the element class of the ghost */
      ghost_class = t8_forest_ghost_get_tree_class (forest, itree);
      ghost_elem = t8_forest_ghost_get_tree_elements (forest, itree);
      num_elements = t8_forest_ghost_tree_num_elements (forest, itree);
      tscheme = t8_forest_get_eclass_scheme (forest, ghost_class);
      for (ielem = 0; ielem < (t8_locidx_t) num_elements; ielem++) {
        elem = t8_element_array_index_locidx (ghost_elem, ielem);
        num_points += tscheme->t8_element_num_corners (elem);
      }
    }
  }
  return num_points;
}

static int
t8_forest_vtk_cells_vertices_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                     t8_tree_t tree,
                                     t8_locidx_t element_index,
                                     t8_element_t *element,
                                     t8_eclass_scheme_c *ts,
                                     int is_ghost,
                                     FILE *vtufile, int *columns,
                                     void **data, T8_VTK_KERNEL_MODUS modus)
{
#if 0
  /* if we eventually implement scaling the elements, activate this line */
  double              midpoint[3];
#endif
  double              element_coordinates[3];
  int                 num_el_vertices, ivertex;
  int                 freturn;
  t8_element_shape_t  element_shape;

  if (modus != T8_VTK_KERNEL_EXECUTE) {
    /* Nothing to do if we are in Init or clean up mode */
    return 1;
  }

  /* TODO: be careful with pyramid class here.
   *       does this work too over tree->class or do we need something else?
   */

#if 0
  /* if we eventually implement scaling the elements, activate this line */
  t8_forest_element_centroid (forest, ltree_id, element, midpoint);
#endif
  element_shape = ts->t8_element_shape (element);
  num_el_vertices = t8_eclass_num_vertices[element_shape];
  for (ivertex = 0; ivertex < num_el_vertices; ivertex++) {
    t8_forest_element_coordinate (forest, ltree_id, element,
                                  t8_eclass_vtk_corner_number[element_shape]
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
t8_forest_vtk_write_vertices (t8_forest_t forest, FILE *vtufile)
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
                                         t8_element_t *elements,
                                         t8_eclass_scheme_c *ts,
                                         int is_ghost,
                                         FILE *vtufile, int *columns,
                                         void **data,
                                         T8_VTK_KERNEL_MODUS modus)
{
  int                 ivertex, num_vertices;
  int                 freturn;
  t8_locidx_t        *count_vertices;
  t8_element_shape_t  element_shape;

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
  element_shape = ts->t8_element_shape (elements);
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
t8_forest_vtk_cells_offset_kernel (t8_forest_t forest, t8_locidx_t ltree_id,
                                   t8_tree_t tree,
                                   t8_locidx_t element_index,
                                   t8_element_t *element,
                                   t8_eclass_scheme_c *ts,
                                   int is_ghost,
                                   FILE *vtufile, int *columns,
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
                                 t8_element_t *element,
                                 t8_eclass_scheme_c *ts,
                                 int is_ghost,
                                 FILE *vtufile, int *columns,
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
                                  t8_element_t *element,
                                  t8_eclass_scheme_c *ts,
                                  int is_ghost,
                                  FILE *vtufile, int *columns,
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
                                 t8_element_t *element,
                                 t8_eclass_scheme_c *ts,
                                 int is_ghost,
                                 FILE *vtufile, int *columns,
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
                                   t8_element_t *element,
                                   t8_eclass_scheme_c *ts,
                                   int is_ghost,
                                   FILE *vtufile, int *columns,
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
                                      t8_element_t *element,
                                      t8_eclass_scheme_c *ts,
                                      int is_ghost,
                                      FILE *vtufile, int *columns,
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
                                   t8_element_t *element,
                                   t8_eclass_scheme_c *ts,
                                   int is_ghost,
                                   FILE *vtufile, int *columns,
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
                                   t8_element_t *element,
                                   t8_eclass_scheme_c *ts,
                                   int is_ghost,
                                   FILE *vtufile, int *columns,
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
                                      t8_element_t *element,
                                      t8_eclass_scheme_c *ts,
                                      int is_ghost,
                                      FILE *vtufile, int *columns,
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
                                      t8_element_t *element,
                                      t8_eclass_scheme_c *ts,
                                      int is_ghost,
                                      FILE *vtufile, int *columns,
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
t8_forest_vtk_write_cell_data (t8_forest_t forest, FILE *vtufile,
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
t8_forest_vtk_write_cells (t8_forest_t forest, FILE *vtufile,
                           int write_treeid,
                           int write_mpirank,
                           int write_level, int write_element_id,
                           int write_ghosts, int num_data,
                           t8_vtk_data_field_t *data)
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
t8_forest_vtk_write_points (t8_forest_t forest, FILE *vtufile,
                            int write_ghosts,
                            int num_data, t8_vtk_data_field_t *data)
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
                          int num_data, t8_vtk_data_field_t *data)
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

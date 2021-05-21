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

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#include <vtkUnsignedCharArray.h>
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default_cxx.hxx>
#include <t8_forest_vtk.h>
#include <t8_vtk.h>
#include <t8_cmesh.h>
#include <t8_element_cxx.hxx>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_vec.h>
#include "t8_cmesh/t8_cmesh_trees.h"


static t8_cmesh_t
t8_step2_build_prismcube_coarse_mesh (sc_MPI_Comm comm)
{
  t8_cmesh_t          cmesh;

  /* Build a coarse mesh of 2 prism trees that form a cube. */
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_PRISM, comm, 0, 0, 0);
  return cmesh;
}

static t8_forest_t
t8_step2_build_uniform_forest (sc_MPI_Comm comm, t8_cmesh_t cmesh, int level)
{
  t8_forest_t         forest;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();
  forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

  return forest;
}

static void
t8_step2_destroy_forest (t8_forest_t forest)
{
  t8_forest_unref (&forest);
}



int main(int argc, char* argv[])
{

  // Parse command line arguments
  /*if (argc != 2)
  {
    std::cout << "Required arguments: OutputFilename.vtu" << std::endl;
    return EXIT_FAILURE;
  }*/

  std::string filename = "write.pvtu";
/*
  Brauchen für das unstructured Grid die Koordinaten der Punkte
  ,die Id der Punkte für die einzelnen Elemente(in vtk Reihenfolge)
  und ein Array der Cell Types also Elementtypen in vtk Cell types 
  konvertiert
  brauchen sowohl loc id als auch glo id der Punkte (konvertierung)
  iteriere für ein forest über alle trees, elemente und punkte und
  greife dann auf die punkte, elementtypen und loc und glo id zu
  sollte dann in eine funktion geschrieben werden, die prefix und
  forest als Input nimmt
  - krieg es für das forest zum laufen
  - wie muss das cellTypeArray gefüllt werden?
   gibt bisher error: invalid user-defined conversion from ‘vtkSmartPointer<vtkUnsignedCharArray>’ to ‘int’
   die funktion ist überladen und es gibt eine version mit int und eine mit vtkUnsignedCharArray
  {VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA} funktioniert auch nicht
  vtkSmartPointer<vtkUnsignedCharArray> cellTypeArray; quad=9, pyramid=14
 */

  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  const int           level = 3;
  t8_locidx_t         local_num_elements;
  t8_gloidx_t         global_num_elements;
  t8_element_t       *ielement;
  t8_tree_t           tree;
  t8_locidx_t         itree, ivertex;
  t8_locidx_t         element_index;
  t8_ctree_t          ctree;
  double             *vertices;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the leg levels. */
  t8_init (SC_LP_PRODUCTION);

 /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /* Create the cmesh from step1 */
  cmesh = t8_step2_build_prismcube_coarse_mesh (comm);

  /* Build the uniform forest, it is automatically partitioned among the processes. */
  forest = t8_step2_build_uniform_forest (comm, cmesh, level);
/*
  cmesh = forest->cmesh;
  int i=0;
  double x=0;
  double y=0;
  double z=0;
  double  coordinates[3];
  double pnt;
  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> cellArray;
for (itree = 0;
    itree < t8_forest_get_num_local_trees (forest); itree++) {
      ctree = t8_cmesh_get_tree (cmesh,
                               t8_forest_ltreeid_to_cmesh_ltreeid (forest,
                                                                   itree));
      vertices = ((double *)
                t8_cmesh_get_attribute (cmesh, t8_get_package_id (), 0,
                                        ctree->treeid));
      tree = t8_forest_get_tree (forest, itree);
  for (ielement = 0;
           ielement < t8_forest_get_tree_num_elements (forest,
                                                       itree);
           ielement++) {
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
        pnt={x,y,z};
        points->InsertPoint(i,x[i]);
        cellArray->InsertNextCell(t8_eclass_num_vertices[tree->eclass],pts[i]);
                                      }

    }
  }*/
  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> cellArray;
  vtkNew<vtkUnsignedCharArray> unsignedCharArray;
  int i;
  static double x[9][3]={{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
  {1,1,0}, {1,0,1}, {1,1,1}, {0,1,1}, {0.5,0.5,2}};
  static vtkIdType pts[7][5]={{0,1,2,4}, {0,1,5,3}, {5,3,7,6},
  {2,4,6,7}, {2,3,7,0}, {1,4,5,6}, {3,5,6,7,8}};
  static const double y[7]={9,9,9,9,9,9,14};
  // Load the point, cell, and data attributes.
  for (i=0; i<10; i++) points->InsertPoint(i,x[i]);
  for (i=0; i<6; i++) cellArray->InsertNextCell(4,pts[i]);

  vtkNew<vtkUnstructuredGrid> unstructuredGrid;
  unstructuredGrid->SetPoints(points);
  for (i=6; i<7; i++) cellArray->InsertNextCell(5,pts[i]);
  for (i=0;i<7;i++) unsignedCharArray->InsertNextTuple(*y[i]);
  unstructuredGrid->SetCells(*unsignedCharArray,cellArray);
  // Write file
  vtkNew<vtkXMLPUnstructuredGridWriter> writer;
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();

  return EXIT_SUCCESS;
}
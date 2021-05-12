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

  std::string filename = "write.vtu";
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
  vtkSmartPointer<vtkUnsignedCharArray> cellTypeArray;
 */

  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> cellArray;
  vtkNew<vtkUnsignedCharArray> cellTypeArray;
  int i;
  static double x[8][3]={{0,0,0}, {2,0,0}, {1,1,0}, {0,1,0},
  {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};
  static vtkIdType pts[6][4]={{0,1,2,3}, {4,5,6,7}, {0,1,5,4},
  {1,2,6,5}, {2,3,7,6}, {3,0,4,7}};
  static vtkCellTypes y[6]={VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA,VTK_TETRA};
  // Load the point, cell, and data attributes.
  for (i=0; i<8; i++) points->InsertPoint(i,x[i]);
  for (i=0; i<6; i++) cellArray->InsertNextCell(4,pts[i]);
  for (i=0; i<6; i++) cellTypeArray->InsertNextTuple(i,y[i]);

  vtkNew<vtkUnstructuredGrid> unstructuredGrid;
  unstructuredGrid->SetPoints(points);
  unstructuredGrid->SetCells(cellTypeArray,cellArray);

  // Write file
  vtkNew<vtkXMLPUnstructuredGridWriter> writer;
  writer->SetFileName(filename.c_str());
  writer->SetInputData(unstructuredGrid);
  writer->Write();

  return EXIT_SUCCESS;
}
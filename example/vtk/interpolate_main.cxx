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

#include <example/vtk/interpolate.hxx>
#include <example/vtk/Callbacks.h>


static void
t8_pipeline (vtkSmartPointer < vtkDataSet >data, sc_MPI_Comm comm)
{

  MeshAdapter *test = new MeshAdapter(comm);

  t8_debugf ("[D] Constructed class\n");

  test->SetLevelMaximum(6);
  test->SetThreshold(0.1);

  /* Set Division Factors */
    test->SetDivisionFactorX(1);
    test->SetDivisionFactorY(1);
    test->SetDivisionFactorZ(1);

    /* Set the element type used for refinement */
  test->SetElementType(1);

  
  test->Initialize(data);

  const char          filename[11] = "class_test";
  const char          first_data[13] = "class_test_0";

  const char         *data_name[3] = {
    "average",
    "data1",
    "data2"
  };
  t8_vtk_data_type_t  types[3];

  types[0] = T8_VTK_SCALAR;
  types[1] = T8_VTK_SCALAR;
  types[2] = T8_VTK_SCALAR;

  test->WritePVTU (first_data, data_name, types);

  for (int i = 0; i < 6; i++) {
    t8_debugf("----------- Round %i -------------\n", i);

    test->Adapt (t8Callbacks::t8_adapt_callback_non_empty, 
                t8Callbacks::t8_itertate_replace_pointids, 0);
    t8_debugf("----------- Adapt %i -------------\n", i);

    test->Partition ();
    t8_debugf("----------- Partition %i -------------\n", i);

    test->SetElements ();
    t8_debugf("----------- SetElements %i -------------\n", i);

    char                fileprefix[12];
    snprintf (fileprefix, BUFSIZ, "%s_%i", filename, i + 1);
    test->WritePVTU (fileprefix, data_name, types);
    t8_debugf("----------- WritePvtu %i -------------\n", i);

  }
  t8_debugf ("[D] wrote vtk file\n");
  return;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  mpiret = sc_MPI_Init (&argc, &argv);
  sc_MPI_Comm         comm = sc_MPI_COMM_WORLD;
  SC_CHECK_MPI (mpiret);
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEBUG);

  vtkSmartPointer < vtkDataSet >vtk_grid = t8_vtk_reader
    ("/group/HPC/Projects/visplore/Examples/GAIA/gaia-parallel_0.pvtu",
     //("/localdata1/knap_da/projects/t8code/t8code/test/testfiles/test_vtk_tri.vtu",
     1, 0, comm, VTK_PARALLEL_UNSTRUCTURED_FILE);

  t8_debugf("[D] datatset has %i points after reading\n", vtk_grid->GetNumberOfPoints());
  t8_debugf ("[D] read successful\n");
  t8_pipeline (vtk_grid, comm);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}

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

#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_vtk_reader.hxx>

#if T8_WITH_VTK
#include <vtkGenericDataObjectReader.h>
#include <vtkPolyData.h>
#endif
T8_EXTERN_C_BEGIN ();

/*Construct a cmesh given a filename and a*/
t8_cmesh_t
t8_cmesh_read_from_vtk (const char *prefix, const int num_files)
{
  t8_cmesh_t          cmesh;
//#if T8_WITH_VTK
  vtkSmartPointer < vtkGenericDataObjectReader > reader =
    vtkSmartPointer < vtkGenericDataObjectReader >::New ();
  reader->SetFileName (prefix);
  reader->Update ();

  if (reader->IsFilePolyData ()) {
    t8_debugf ("[D] reader got polyData\n");
    vtkPolyData        *output = reader->GetPolyDataOutput ();
    t8_debugf ("[D] reader has read %lli points",
               output->GetNumberOfPoints ());
  }
  else {
    t8_global_errorf ("Data in file is not polydata");
  }
//#else
  /*TODO: Proper return value to prevent compiler-errors */
  t8_global_errorf
    ("WARNING: t8code is not linked againt the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
//#endif
  return cmesh;
}

T8_EXTERN_C_END ();
